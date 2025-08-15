classdef StatisticalAnalyzer < handle
    properties (Access = private)
        results_data cell
        results_table table
    end
    
    methods
        function obj = StatisticalAnalyzer(varargin)            
            p = inputParser;
            addOptional(p, 'results', {}, @(x) iscell(x));
            addParameter(p, 'results_file', '', @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});
            
            results = p.Results.results;
            results_file = p.Results.results_file;
            
            if ~isempty(results) && ~isempty(results_file)
                error('Cannot provide both results and results_file');
            elseif isempty(results) && isempty(results_file)
                error('Provide either results or results_file');
            elseif ~isempty(results)
                obj.results_data = results;
                obj.results_table = obj.convert_results_to_table(results);
            else
                loaded_data = load(results_file);
                field_names = fieldnames(loaded_data);
                obj.results_data = loaded_data.(field_names{1});
                obj.results_table = obj.convert_results_to_table(obj.results_data);
            end
        end
        
        function results_table = convert_results_to_table(~, results)
            n_results = length(results);
            
            C_values = zeros(n_results, 1);
            K_values = zeros(n_results, 1);
            simulation_ids = zeros(n_results, 1);
            num_simulations = zeros(n_results, 1);
            
            average_num_users_cell = cell(n_results, 1);
            loss_probabilities_cell = cell(n_results, 1);
            num_arrivals_cell = cell(n_results, 1);
            num_losses_cell = cell(n_results, 1);
            state_probabilities_cell = cell(n_results, 1);
            state_residence_times_cell = cell(n_results, 1);
            
            for i = 1:n_results
                result = results{i};
                C_values(i) = result.config.C;
                K_values(i) = result.config.K;
                simulation_ids(i) = result.config.simulation_id;
                num_simulations(i) = result.config.num_simulations;
                
                average_num_users_cell{i} = result.average_num_users;
                loss_probabilities_cell{i} = result.loss_probabilities;
                num_arrivals_cell{i} = result.num_arrivals;
                num_losses_cell{i} = result.num_losses;
                state_probabilities_cell{i} = result.state_probabilities;
                state_residence_times_cell{i} = result.state_residence_times;
            end
            
            results_table = table(C_values, K_values, simulation_ids, num_simulations, ...
                average_num_users_cell, loss_probabilities_cell, ...
                num_arrivals_cell, num_losses_cell, ...
                state_probabilities_cell, state_residence_times_cell, ...
                'VariableNames', {'C', 'K', 'simulation_id', 'num_simulations', ...
                'average_num_users', 'loss_probabilities', ...
                'num_arrivals', 'num_losses', ...
                'state_probabilities', 'state_residence_times'});
        end
        
        function mean_table = calculate_mean_and_group_by_ck(obj, column_name, alias)            
            if ~ismember(column_name, obj.results_table.Properties.VariableNames)
                error('Column %s not found in results table', column_name);
            end
            exploded_table = obj.explode_column(obj.results_table, column_name);            
            [unique_ck, ~, group_idx] = unique(exploded_table(:, {'C', 'K'}));            
            mean_values = splitapply(@mean, exploded_table.(column_name), group_idx);            
            mean_table = table(unique_ck.C, unique_ck.K, mean_values, ...
                'VariableNames', {'C', 'K', alias});
        end
        
        function ci_table = calculate_confidence_intervals_df(obj, c_value, confidence_level)            
            if nargin < 3
                confidence_level = 0.95;
            end            
            c_filtered = obj.results_table(obj.results_table.C == c_value, :);
            
            if isempty(c_filtered)
                error('No data found for C = %d', c_value);
            end            
            exploded_table = obj.explode_column(c_filtered, 'average_num_users');            
            [unique_k, ~, group_idx] = unique(exploded_table.K);            
            means = splitapply(@mean, exploded_table.average_num_users, group_idx);
            stds = splitapply(@(x) std(x, 1), exploded_table.average_num_users, group_idx); % ddof=1
            ns = splitapply(@length, exploded_table.average_num_users, group_idx);
            
            alpha = (1 - confidence_level) / 2;
            df = ns - 1;
            t_critical = tinv(1 - alpha, df);
            margin_error = t_critical .* (stds ./ sqrt(ns));
            
            ci_lower = means - margin_error;
            ci_upper = means + margin_error;
            
            confidence_levels = repmat(confidence_level, length(unique_k), 1);
            C_values = repmat(c_value, length(unique_k), 1);
            
            ci_table = table(unique_k, means, stds, ns, ci_lower, ci_upper, ...
                confidence_levels, C_values, ...
                'VariableNames', {'K', 'mean', 'std', 'n', 'ci_lower', 'ci_upper', ...
                'confidence_level', 'C'});
        end   

        function exploded_table = explode_column(~, input_table, column_name)
            n_rows = height(input_table);
            all_values = [];
            all_k = [];
            all_c = [];
            for i = 1:n_rows
                values = input_table.(column_name){i};
                n_values = length(values);
                all_values = [all_values; values(:)];
                all_k = [all_k; repmat(input_table.K(i), n_values, 1)];
                all_c = [all_c; repmat(input_table.C(i), n_values, 1)];
            end
            exploded_table = table(all_c, all_k, all_values, ...
                'VariableNames', {'C', 'K', column_name});
        end

        function plot_mean_grouped_by_ck(~, mean_table, alias)
            fig_name = sprintf('Mean %s by C and K', strrep(alias, '_', ' '));
            figure('Name', fig_name, 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
            clf;
            hold on;
            
            unique_c = unique(mean_table.C);
            colors = lines(length(unique_c));
            
            for i = 1:length(unique_c)
                c_val = unique_c(i);
                c_mask = mean_table.C == c_val;
                c_data = mean_table(c_mask, :);                
                c_data = sortrows(c_data, 'K');
                
                plot(c_data.K, c_data.(alias), '-o', ...
                    'LineWidth', 2, 'MarkerSize', 8, ...
                    'Color', colors(i, :), ...
                    'DisplayName', sprintf('C=%d', c_val));
            end
            
            xlabel('Buffer Length (K)', 'FontSize', 12);
            ylabel(strrep(alias, '_', ' '), 'FontSize', 12);
            title(strrep(alias, '_', ' '), 'FontSize', 14);
            grid on;
            legend('show', 'Location', 'best');
            hold off;
        end

        function plot_confidence_intervals(~, ci_table, print_statistics)
            if nargin < 3
                print_statistics = true;
            end
            
            K_values = ci_table.K;
            means = ci_table.mean;
            ci_lower = ci_table.ci_lower;
            ci_upper = ci_table.ci_upper;
            
            lower_errors = means - ci_lower;
            upper_errors = ci_upper - means;
            
            figure('Position', [100, 100, 900, 700]);
            errorbar(K_values, means, lower_errors, upper_errors, ...
                '-', 'LineWidth', 2, 'MarkerSize', 8, ...
                'CapSize', 10);
            
            xlabel('K', 'FontSize', 12);
            ylabel('Average number of users', 'FontSize', 12);
            title(sprintf('Average number of users (CI) - C = %d', ci_table.C(1)), 'FontSize', 14);
            grid on;
            xticks(K_values);            
            y_min = floor(min(ci_lower));
            y_max = ceil(max(ci_upper));
            yticks(y_min:2:y_max);
            
            if print_statistics
                fprintf('\nConfidence Interval Statistics for C = %d:\n', ci_table.C(1));
                fprintf('Confidence Level: %.1f%%\n', ci_table.confidence_level(1) * 100);
                fprintf('\nK\tMean\t\tCI Lower\tCI Upper\tWidth\t\tStd Dev\t\tn\n');
                fprintf('%s\n', repmat('-', 1, 80));
                
                for i = 1:height(ci_table)
                    width = ci_table.ci_upper(i) - ci_table.ci_lower(i);
                    fprintf('%d\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%d\n', ...
                        ci_table.K(i), ci_table.mean(i), ci_table.ci_lower(i), ...
                        ci_table.ci_upper(i), width, ci_table.std(i), ci_table.n(i));
                end
            end
        end
    end
end