function run_single()
    clear; %#ok<CLEAR0ARGS>
    clc; close all;    
    lambda = 1/10;
    rho = 1.05;
    base_config = SimulationConfiguration( ...
        1000, ...             % simulation_id
        lambda, ...          % lambda_
        rho, ...                     % rho
        2, ...                       % C
        30, ...                      % K
        100, ...                     % num_simulations
        50000, ...                 % length_simulation
        1000000);                       % max_iterations
    simulator = E2E2CKSimulator(base_config);
    fprintf("Running a simulation..\n")
    tic;
    simulator.simulate();
    duration = toc;
    fprintf('Simulation completed in %.2f seconds.\n', duration);
    
    % result = simulator.result;
    % fprintf('\nSimulation Results Summary:\n');
    % fprintf('Configuration: C=%d, K=%d, lambda=%.3f, rho=%.3f\n', ...
    %     result.config.C, result.config.K, result.config.lambda, result.config.rho);
    % fprintf('Number of simulations: %d\n', length(result.num_arrivals));
    % fprintf('Average arrivals per simulation: %.2f\n', mean(result.num_arrivals));
    % fprintf('Average losses per simulation: %.2f\n', mean(result.num_losses));
    % fprintf('Average loss probability: %.4f\n', mean(result.loss_probabilities));
    % fprintf('Average number of users: %.4f\n', mean(result.average_num_users));
end

function configs = create_config_grid_ck(base_config, c_values, k_max)
    total_configs = sum(arrayfun(@(c) (k_max - c + 1), c_values));

    configs = cell(1, total_configs);
    simulation_id = base_config.simulation_id;
    index = 1;

    for i = 1:length(c_values)
        c = c_values(i);
        k_values = c:k_max;

        for j = 1:length(k_values)
            k = k_values(j);

            config = SimulationConfiguration( ...
                simulation_id, ...
                base_config.lambda_arrival, ...
                base_config.rho, ...
                base_config.num_arrival_stages, ...
                base_config.num_service_stages, ...
                c, ...
                k, ...
                base_config.num_simulations, ...
                base_config.length_simulation);

            configs{index} = config;
            index = index + 1;
            simulation_id = simulation_id + 1;
        end
    end
end


function results = run_simulations(configs)    
    results = cell(length(configs), 1);
    
    fprintf('Running %d simulations...\n', length(configs));
    
    for i = 1:length(configs)
        try
            config = configs{i};
            fprintf('Starting simulation: C=%d, K=%d, ID=%d\n', ...
                config.C, config.K, config.simulation_id);
            
            simulator = VectorizedSimulator(config);
            simulator.simulate();
            results{i} = simulator.result;
            
            fprintf('Completed simulation: C=%d, K=%d, ID=%d\n', ...
                config.C, config.K, config.simulation_id);
        catch ME
            fprintf('Error in simulation C=%d, K=%d, ID=%d: %s\n', ...
                config.C, config.K, config.simulation_id, ME.message);
            results{i} = struct('error', ME.message);
        end
    end
end

run_single();