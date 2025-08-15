function run_single() %#ok<DEFNU>
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
                base_config.lambda, ...
                base_config.rho, ...
                c, ...
                k, ...
                base_config.num_simulations, ...
                base_config.length_simulation, ...
                base_config.max_iterations);

            configs{index} = config;
            index = index + 1;
            simulation_id = simulation_id + 1;
        end
    end
end

function run_parallel()
    clear; %#ok<CLEAR0ARGS>
    clc; close all;
    
    % Ensure we have exactly 5 workers
    current_pool = gcp('nocreate');
    if isempty(current_pool)
        parpool('local', 5);
        fprintf('Started parallel pool with 5 workers.\n');
    elseif current_pool.NumWorkers ~= 5
        delete(current_pool);
        parpool('local', 5);
        fprintf('Restarted parallel pool with 5 workers.\n');
    else
        fprintf('Using existing parallel pool with %d workers.\n', current_pool.NumWorkers);
    end
    
    % Base configuration
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
        1000000);                    % max_iterations
    
    c_values = [2, 3];
    k_max = 30;
    configs = create_config_grid_ck(base_config, c_values, k_max);
    fprintf('Created %d configurations to simulate.\n', length(configs));
    fprintf('Running parallel simulations...\n');
    results = cell(size(configs));
    tic;
    parfor i = 1:length(configs)
        config = configs{i};
        simulator = E2E2CKSimulator(config);
        simulator.simulate();
        results{i} = simulator.result;
        
        fprintf('Worker completed simulation %d: C=%d, K=%d\n', ...
            config.simulation_id, config.C, config.K);
    end
    duration = toc;
    
    fprintf('\nAll %d simulations completed in %.2f seconds.\n', length(configs), duration);
    % Display summary results
    fprintf('\nSummary of Results:\n');
    fprintf('%-5s %-5s %-10s %-12s %-15s %-15s %-15s\n', 'C', 'K', 'Arrivals', 'Losses', 'Loss Prob', 'Avg Users', 'Total Iterations');
    fprintf('%-5s %-5s %-10s %-12s %-15s %-15s %-15s\n', '-', '-', '--------', '------', '---------', '---------', '-----------------');

    for i = 1:length(results)
        result = results{i};
        fprintf('%-5d %-5d %-10.0f %-12.0f %-15.4f %-15.4f %-15d\n', ...
            result.config.C, result.config.K, ...
            mean(result.num_arrivals), mean(result.num_losses), ...
            mean(result.loss_probabilities), mean(result.average_num_users), ...
            mean(result.total_iterations));
    end
end

run_parallel();