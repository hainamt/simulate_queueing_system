function main()
    clear; %#ok<CLEAR0ARGS>
    clc; close all;    
    lambda = 1/10;
    rho = 1.05;
    C = 2;
    K = 30;
    num_simulations = 100;
    length_simulation = 50000;
    max_iterations = 1000000;
    num_workers = 1;
    base_config = SimulationConfiguration( ...
        1000, ...             % simulation_id
        lambda, ...
        rho, ...
        C, ...
        K, ...
        num_simulations, ...
        length_simulation, ...
        max_iterations);

    configs = create_config_grid_ck(base_config, [2, 3], 30);
    results = run_parallel(configs, num_workers);

    analyzer = StatisticalAnalyzer('results', results);
    mean_num_user_table = analyzer.calculate_mean_and_group_by_ck('average_num_users', 'mean_average_num_users');
    mean_loss_prob_table = analyzer.calculate_mean_and_group_by_ck('loss_probabilities', 'mean_loss_probabilities');
    ci_table = analyzer.calculate_confidence_intervals_df(2, 0.95);

    analyzer.plot_mean_grouped_by_ck(mean_num_user_table, 'mean_average_num_users');
    analyzer.plot_mean_grouped_by_ck(mean_loss_prob_table, 'mean_loss_probabilities');
    analyzer.plot_confidence_intervals(ci_table, false);
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

function results = run_parallel(configs, num_workers)
    arguments
        configs (1,:) cell
        num_workers uint16
    end
    clc; close all;

    current_pool = gcp('nocreate');
    if isempty(current_pool)
        parpool('local', num_workers);
        fprintf('Started parallel pool with %d workers.\n', num_workers);
    elseif current_pool.NumWorkers ~= num_workers
        delete(current_pool);
        parpool('local', num_workers);
        fprintf('Restarted parallel pool with %d workers.\n', num_workers);
    else
        fprintf('Using existing parallel pool with %d workers.\n', current_pool.NumWorkers);
    end
    
    fprintf('Simulate %d configurations.\n', length(configs));
    results = cell(size(configs));
    tic;
    parfor i = 1:length(configs)
        config = configs{i};
        simulator = E2E2CKSimulator(config);
        simulator.simulate();
        results{i} = simulator.result;

        fprintf('Completed simulation %d: C=%d, K=%d\n', ...
            config.simulation_id, config.C, config.K);
    end
    duration = toc;
    fprintf('\nAll %d simulations completed in %.2f seconds.\n', length(configs), duration);
end


%    fprintf('\nSummary of Results:\n');
%    fprintf('%-5s %-5s %-10s %-12s %-15s %-15s %-15s\n', 'C', 'K', 'Arrivals', 'Losses', 'Loss Prob', 'Avg Users', 'Total Iterations');
%    fprintf('%-5s %-5s %-10s %-12s %-15s %-15s %-15s\n', '-', '-', '--------', '------', '---------', '---------', '-----------------');
%
%    for i = 1:length(results)
%        result = results{i};
%        fprintf('%-5d %-5d %-10.0f %-12.0f %-15.4f %-15.4f %-15d\n', ...
%            result.config.C, result.config.K, ...
%            mean(result.num_arrivals), mean(result.num_losses), ...
%            mean(result.loss_probabilities), mean(result.average_num_users), ...
%            mean(result.total_iterations));
%    end