classdef VectorizedE2E2CKSimulator < handle
    properties
        config SimulationConfiguration
        C uint8 = 2
        K uint8 = 2

        lambda double
        mu double
        rho double

        num_simulations uint16
        length_simulation double

        num_arrivals (:,1) uint32
        num_losses (:,1) uint32
        x_in (:,1) double
        num_event_types (1,1) uint32
        event_table (:,:) double
        number_users_when_event (:,:) uint32
        clk (:,1) double
        previous_queue_change_time (:,1) double

        total_customers (:,1) uint32
        server_states (:,:) uint8
        state_residence_time (:,:) double
        actual_end_times (:,1) double
        active_simulations (:,1) logical
        rng

        result VectorizedSimulationResult
    end

    methods
        function obj = VectorizedE2E2CKSimulator(config)
            arguments
                config SimulationConfiguration
            end
            obj.config = config;
            obj.C = config.C;
            obj.K = config.K;

            obj.lambda = config.lambda;
            obj.mu = config.mu;
            obj.rho = config.rho;

            obj.num_simulations = config.num_simulations;
            obj.length_simulation = config.length_simulation;

            obj.num_arrivals = zeros(obj.num_simulations, 1);
            obj.num_losses = zeros(obj.num_simulations, 1);

            obj.x_in = ones(obj.num_simulations, 1);
            obj.num_event_types = 2 + 2 * obj.C;
            obj.event_table = inf([obj.num_simulations, obj.num_event_types]);

            obj.clk = zeros(obj.num_simulations, 1);
            obj.previous_queue_change_time = zeros(obj.num_simulations, 1);

            obj.total_customers = zeros(obj.num_simulations, 1);
            obj.server_states = zeros(obj.num_simulations, obj.C);
            obj.state_residence_time = zeros(obj.num_simulations, obj.K + 1);
            obj.actual_end_times = repmat(obj.length_simulation, obj.num_simulations, 1);

            obj.active_simulations = true(obj.num_simulations, 1);
            
            obj.rng = RandStream('mt19937ar', 'Seed', round(posixtime(datetime("now"))));

        end

        %% utils
        function exponential_samples = rand_exp(obj, sim_indices, rate)
            size = length(sim_indices);
            uniform_samples = rand(obj.rng, size, 1);
            uniform_samples(uniform_samples == 0) = rand(obj.rng, sum(uniform_samples == 0), 1);
            exponential_samples = -log(1 - uniform_samples) / rate;
            
            if size == 1
                exponential_samples = exponential_samples(1);
            end
        end

        function update_state_residence_time(obj, sim_indices)
            time_diff = obj.clk(sim_indices) - obj.previous_queue_change_time(sim_indices);
            current_states = obj.total_customers(sim_indices) + 1;
            obj.state_residence_time(sim_indices, current_states) = obj.state_residence_time(sim_indices, current_states) + time_diff;
        end
        
        %% handler
        function handle_a1_event(obj, sim_indices)
            obj.x_in(sim_indices) = 2;
            obj.event_table(sim_indices, 1) = inf;
            obj.event_table(sim_indices, 2) = obj.clk(sim_indices) + obj.rand_exp(sim_indices, 2*obj.lambda);
        end

        function handle_a2_event(obj, sim_indices)
            obj.num_arrivals(sim_indices) = obj.num_arrivals(sim_indices) + 1;
            can_enter = obj.total_customers(sim_indices) < obj.K;
            entering_sims = sim_indices(can_enter);

            if ~isempty(entering_sims)
                obj.update_state_residence_time(entering_sims);
                obj.total_customers(entering_sims) = obj.total_customers(entering_sims) + 1;

                server_states_entering = obj.server_states(entering_sims);
                idle_mask = server_states_entering == 0;
                has_idle_server = any(idle_mask, 2);
                server_assignments = -ones([length(entering_sims), 1]);

                if any(has_idle_server)
                    [~, first_idle_indices] = max(idle_mask, [], 2);
                    server_assignments(has_idle_server) = first_idle_indices(has_idle_server);

                    valid_mask = server_assignments >= 0;
                    valid_entering_sims = entering_sims(valid_mask);
                    valid_server_assignments = server_assignments(valid_mask);

                    if ~isempty(valid_server_assignments)
                        obj.server_states(valid_entering_sims, valid_server_assignments) = 1;
                        s1_event_indices = 2 + 2 * valid_server_assignments - 1;
                        service_times = obj.rand_exp(valid_entering_sims, 2* obj.mu);
                        obj.event_table(valid_entering_sims, s1_event_indices) = obj.clk(valid_entering_sims) + service_times;
                    end
                end
                obj.previous_queue_change_time(entering_sims) = obj.clk(entering_sims);
            end

            blocked = obj.total_customers(sim_indices) >= obj.K;
            blocked_sims = sim_indices(blocked);
            if ~isempty(blocked_sims)
                obj.num_losses(blocked_sims) = obj.num_losses(blocked_sims) + 1;
            end
            obj.x_in(sim_indices) = 1;
            obj.event_table(sim_indices, 1) = obj.clk(sim_indices) + obj.rand_exp(sim_indices, 2*obj.lambda);
            obj.event_table(sim_indices, 2) = inf;
        end

        function handle_s1_event(obj, sim_indices, server_id)
            obj.server_states(sim_indices, server_id) = 2;
            s1_event_indices = 2 + 2 * server_ids - 1;
            s2_event_indices = s1_event_indices + 1;

            obj.event_table(sim_indices, s1_event_indices) = inf;
            service_times = obj.rand_exp(sim_indices, 2* obj.mu);
            obj.event_table(sim_indices, s2_event_indices) = obj.clk(sim_indices) + service_times;
        end

        function handle_s2_event(obj, sim_indices, server_id)
            obj.update_state_residence_time(sim_indices);
            obj.total_customers(sim_indices) = obj.total_customers(sim_indices) - 1;

            s1_event_indices = 2 + 2 * server_id - 1;
            s2_event_indices = s1_event_indices + 1;

            becomes_idle = obj.total_customers(sim_indices) < obj.C;
            continues_serving = ~becomes_idle;

            idle_sims = sim_indices(becomes_idle);
            idle_servers = server_ids(becomes_idle);
            if ~isempty(idle_sims)
                obj.server_states(idle_sims, idle_servers) = 0;
                idle_s1_indices = 2 + 2 * idle_servers - 1;
                obj.event_table(idle_sims, idle_s1_indices) = inf;
            end

            serving_sims = sim_indices(continues_serving);
            serving_servers = server_ids(continues_serving);
            if ~isempty(serving_sims)
                obj.server_states(serving_sims, serving_servers) = 1;
                serving_s1_indices = 2 + 2 * serving_servers - 1;
                service_times = obj.rand_exp(serving_sims, 2 * obj.mu);
                obj.event_table(serving_sims, serving_s1_indices) = obj.clk(serving_sims) + service_times;
            end

            obj.event_table(sim_indices, s2_event_indices) = inf;
            obj.previous_queue_change_time(sim_indices) = obj.clk(sim_indices);
        end

        %% simulation
        function simulate(obj)
            % setup first event
            obj.event_table(obj.active_simulations, 1) = obj.clk(obj.active_simulations) + obj.rand_exp(obj.active_simulations, 2 * obj.lambda);
            fprintf('Initial event table:\n');
            disp(obj.event_table);
            iteration = 0;
            max_iterations = 100;

            while any(obj.active_simulations) && iteration < max_iterations
                % current iteration
                fprintf('Current iteration: %d\n', iteration);

                indices_current_simulations = find(obj.active_simulations);
                fprintf('All Active Simulations indices: %s\n', mat2str(indices_current_simulations));

                if isempty(indices_current_simulations)
                    break;
                end
                
                [next_event_times, next_event_indices] = min(obj.event_table(indices_current_simulations, :), [], 2);
                obj.clk(indices_current_simulations) = next_event_times;

                is_finished = obj.clk(indices_current_simulations) >= obj.length_simulation;
                indices_finished_sims = indices_current_simulations(is_finished);

                % check if there is finished simulations
                if ~isempty(indices_finished_sims)
                    final_time_diff = obj.actual_end_times(indices_finished_sims) - obj.previous_queue_change_time(indices_finished_sims);
                    final_states = obj.total_customers(indices_finished_sims) + 1;
                    linear_indices = sub2ind(size(obj.state_residence_time), indices_finished_sims, final_states);
                    obj.state_residence_time(linear_indices) = obj.state_residence_time(linear_indices) + final_time_diff;
                end

                obj.active_simulations(indices_finished_sims) = false;
                still_active = indices_current_simulations(~is_finished);
                if isempty(still_active)
                    break;
                end

                still_active_event_indices = next_event_indices(~is_finished);

                for event_type = 1:obj.num_event_types
                    event_mask = still_active_event_indices == event_type;
                    if ~any(event_mask)
                        continue;
                    end

                    % events that are in current event_type
                    event_sims = still_active(event_mask);
                    fprintf('Processing event type %d for simulations: %s\n', event_type, mat2str(event_sims));
                    if event_type == 1
                        obj.handle_a1_event(event_sims);
                        fprintf('Event table after handling A1 events:\n');
                        disp(obj.event_table);
                    elseif event_type == 2
                        obj.handle_a2_event(event_sims);
                        fprintf('Event table after handling A2 events:\n');
                        disp(obj.event_table);
                    else
                        server_of_event = floor((event_type - 3) / 2) + 1;
                        if mod(event_type, 2) == 1
                            obj.handle_s1_event(event_sims, server_of_event);
                            fprintf('Event table after handling S1 events:\n');
                            disp(obj.event_table);
                        else
                            obj.handle_s2_event(event_sims, server_of_event);
                            fprintf('Event table after handling S2 events:\n');
                            disp(obj.event_table);
                        end
                    end
                end

                iteration = iteration + 1;
            end

            obj.generate_statistics();
            fprintf('Finished in %d iterations\n', iteration);
        end

        function generate_statistics(obj)
            fprintf('Generating statistics...\n');
            valid_arrivals = obj.num_arrivals > 0;
            loss_probabilities = zeros(obj.num_simulations, 1);
            loss_probabilities(valid_arrivals) = obj.num_losses(valid_arrivals) ./ obj.num_arrivals(valid_arrivals);
        
            state_probabilities = obj.state_residence_time ./ sum(obj.state_residence_time, 2);
            states = double(0:obj.K)';
            average_num_users = sum(state_probabilities .* states', 2);
        
            obj.result = VectorizedSimulationResult(obj.config, ...
                obj.num_arrivals, ...
                obj.num_losses, ...
                average_num_users, ...
                loss_probabilities, ...
                obj.state_residence_time, ...
                state_probabilities, ...
                repmat(obj.length_simulation, obj.num_simulations, 1));
        end
    end
end