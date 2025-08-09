classdef VectorizedSimulator
    properties
        Property1
        config SimulationConfiguration
        C uint8 = 2
        K uint8 = 2

        lambda double
        mu double
        rho double

        num_arrival_stages uint8
        num_service_stages uint8

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
        rng_seeds (:,1) double
        rng (1,1) RandStream

        result VectorizedSimulationResult
    end

    methods
        function obj = VectorizedSimulator(config)
            arguments
                config SimulationConfiguration
            end
            obj.config = config;
            obj.C = config.C;
            obj.K = config.K;

            obj.lambda = config.lambda;
            obj.mu = config.mu;
            obj.rho = config.rho;

            obj.num_arrival_stages = config.num_arrival_stages;
            obj.num_service_stages = config.num_service_stages;

            obj.num_simulations = config.num_simulations;
            obj.length_simulation = config.length_simulation;

            obj.num_arrivals = zeros(obj.num_simulations, 1);
            obj.num_losses = zeros(obj.num_simulations, 1);

            obj.x_in = ones(obj.num_simulations, 1);
            obj.num_event_types = obj.num_arrival_stages + obj.num_service_stages * obj.C;
            obj.event_table = inf([obj.num_simulations, obj.num_event_types]);
            obj.number_users_when_event = zeros(obj.num_simulations, obj.length_simulation);

            obj.clk = zeros(obj.num_simulations, 1);

            obj.total_customers = zeros(obj.num_simulations, 1);
            obj.server_states = zeros(obj.num_simulations, obj.C);
            obj.state_residence_time = zeros(obj.num_simulations, obj.K + 1);
            obj.actual_end_times = repmat(obj.length_simulation, obj.num_simulations, 1);

            obj.active_simulations = true(obj.num_simulations, 1);
            obj.rng_seeds = randi(1e6, obj.num_simulations, 1);
            obj.rng = RandStream('mt19937ar', 'Seed', sum(100*clock));

        end

        function exponential_samples = rand_exp(sim_indices, rate)
            % Validate input
            if rate <= 0
                error('Rate must be positive.');
            end
        
            size = length(sim_indices);
            uniform_samples = rand(size, 1);       
            uniform_samples(uniform_samples == 0) = rand(sum(uniform_samples == 0), 1);        
            exponential_samples = -log(1 - uniform_samples) / rate;
        
            if size == 1
                exponential_samples = exponential_samples(1);
            end
        end

        function update_state_residence_time(obj, sim_indices)
            time_diff = obj.clk(sim_indices);
            current_states = obj.total_customers(sim_indices);
            obj.state_residence_time(sim_indices, current_states) = obj.state_residence_time(sim_indices, current_states) + time_diff;
        end

        function handle_a1_event(obj, sim_indices)
            obj.x_in(sim_indices) = 2;
            obj.event_table(sim_indices, 1) = inf;
            obj.event_table(sim_indices, 2) = obj.clk(sim_indices) + obj.rand_exp(sim_indices, 2*obj.lambda);
        end

        function handle_a2_event(obj, sim_indices)
            obj.num_arrivals(sim_indices) = obj.num_arrivals(sim_indices) + 1;
            can_enter = obj.total_customers(sim_indices) < obj.K;
            entering_sims = sim_indices(can_enter);

            if size(entering_sims) > 0
                obj.update_state_residence_time(entering_sims);
                obj.total_customers(entering_sims) = obj.total_customers(entering_sims) + 1;

                server_states_entering = obj.server_states(entering_sims);
                idle_mask = server_states_entering == 0;
                has_idle_server = any(idle_mask, 2);
                server_assignments = -ones([size(entering_sims), 1]);

                if any(has_idle_server)
                    [~, first_idle_indices] = max(idle_mask, [], 2);
                    server_assignments(has_idle_server) = first_idle_indices(has_idle_server);

                    valid_mask = server_assignments >= 0;
                    valid_entering_sims = entering_sims(valid_mask);
                    valid_server_assignments = server_assignments(valid_mask);

                    if size(valid_server_assignments) > 0
                        obj.server_states(valid_entering_sims, valid_server_assignments) = 1;
                        s1_event_indices = obj.num_arrival_stages + 2 * valid_server_assignments;
                        service_times = obj.rand_exp(valid_entering_sims, 2* obj.mu);
                        obj.event_table(valid_entering_sims, s1_event_indices) = obj.clk(valid_entering_sims) + service_times;
                    end
                end
                obj.previous_queue_change_time(entering_sims) = obj.clk(entering_sims);
            end

            blocked = obj.total_customers(sim_indices) >= obj.K;
            blocked_sims = sim_indices(blocked);
            if size(blocked_sims) > 0
                obj.num_losses = obj.num_losses + 1;
            end
            obj.x_in(sim_indices) = 1;
            obj.event_table(sim_indices, 1) = obj.clk(sim_indices) + obj.rand_exp(sim_indices, 2*obj.lambda);
            obj.event_table(sim_indices, 2) = inf;
        end

        function handle_s1_event(obj, sim_indices, server_ids)
            obj.server_states(sim_indices, server_ids) = 2;
            s1_event_indices = 2 + 2 * server_ids;
            s2_event_indices = s1_event_indices + 1;

            obj.event_table(sim_indices, s1_event_indices) = inf;
            service_times = obj.rand_exp(sim_indices, 2* obj.mu);
            obj.event_table(valid_entering_sims, s2_event_indices) = obj.clk(valid_entering_sims) + service_times;
        end

        function handle_s2_event(obj, sim_indices, server_ids)
            obj.update_state_residence_time(sim_indices);
            obj.total_customers(sim_indices) = obj.total_customers(sim_indices) - 1;

            s1_event_indices = 2 + 2 * server_ids;
            s2_event_indices = s1_event_indices + 1;

            becomes_idle = obj.total_customers(sim_indices) < obj.C;
            continues_serving = ~becomes_idle;

            idle_sims = sim_indices(becomes_idle);
            idle_servers = server_ids(becomes_idle);
            if size(idle_sims) > 0
                obj.server_states(idle_sims, idle_servers) = 0;
                idle_s1_indices = 2 + 2 * idle_servers;
                obj.event_table(idle_sims, idle_s1_indices) = inf;
            end

            serving_servers = server_ids(continues_serving);
            if size(serving_sims) > 0
                obj.server_states(serving_sims, serving_servers) = 1;
                serving_s1_indices = 2 + 2 * serving_servers;
                service_times = obj.rand_exp(serving_sims, 2 * obj.mu);
                obj.event_table(serving_sims, serving_s1_indices) = obj.clk(serving_sims) + service_times;
            end

            obj.event_table(sim_indices, s2_event_indices) = inf;
            obj.previous_queue_change_time(sim_indices) = obj.clk(sim_indices);
        end

        function simulate(obj)
            active_sims = (1:obj.num_simulations)';
            obj.event_table(active_sims, 1) = obj.clk(active_sims) + obj.rand_exp(active_sims, 2 * obj.lambda);

            iteration = 0;
            max_iterations = 1000000;

            while any(obj.active_simulations) && iteration < max_iterations
                active_sims = find(obj.active_simulations);
                if isempty(active_sims)
                    break;
                end

                [next_event_times, next_event_indices] = min(obj.event_table(active_sims, :), [], 2);
                obj.clk(active_sims) = next_event_times;

                finished = obj.clk(active_sims) >= obj.length_simulation;
                finished_sims = active_sims(finished);

                if ~isempty(finished_sims)
                    final_time_diff = obj.actual_end_times(finished_sims) - obj.previous_queue_change_time(finished_sims);
                    final_states = obj.total_customers(finished_sims) + 1; % MATLAB indexing
                    linear_indices = sub2ind(size(obj.state_residence_time), finished_sims, final_states);
                    obj.state_residence_time(linear_indices) = obj.state_residence_time(linear_indices) + final_time_diff;
                end

                obj.active_simulations(finished_sims) = false;
                still_active = active_sims(~finished);
                if isempty(still_active)
                    break;
                end

                still_active_event_indices = next_event_indices(~finished);

                for event_type = 1:obj.num_event_types
                    event_mask = still_active_event_indices == event_type;
                    if ~any(event_mask)
                        continue;
                    end

                    event_sims = still_active(event_mask);
                    if event_type == 1
                        obj.handle_a1_event(event_sims);
                    elseif event_type == 2
                        obj.handle_a2_event(event_sims);
                    else
                        server_id = ceil((event_type - 2) / 2);
                        if mod(event_type, 2) == 1 % odd event_type (s1 events)
                            server_ids = repmat(server_id, length(event_sims), 1);
                            obj.handle_s1_event(event_sims, server_ids);
                        else % even event_type (s2 events)
                            server_ids = repmat(server_id, length(event_sims), 1);
                            obj.handle_s2_event(event_sims, server_ids);
                        end
                    end
                end

                iteration = iteration + 1;
            end

            obj.generate_statistics();
            fprintf('Finished in %d iterations\n', iteration);
        end

        function obj = generate_statistics(obj)
            valid_arrivals = obj.num_arrivals > 0;
            loss_probabilities = zeros(obj.num_simulations, 1);
            loss_probabilities(valid_arrivals) = obj.num_losses(valid_arrivals) ./ obj.num_arrivals(valid_arrivals);
        
            state_probabilities = obj.state_residence_time ./ sum(obj.state_residence_time, 2);
            states = (0:obj.K)'; % Column vector from 0 to K
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