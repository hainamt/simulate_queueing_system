classdef StatisticalAnalyzer < handle    
    properties
        confidence_level double
        results_table table
        results_cell cell
    end
    
    methods
        function obj = StatisticalAnalyzer(confidence_level, varargin)            
            obj.confidence_level = confidence_level;
            
            p = inputParser;
            addParameter(p, 'results', []);
            addParameter(p, 'results_file', '');
            parse(p, varargin{:});
            
            results = p.Results.results;
            results_file = p.Results.results_file;
            
            if ~isempty(results) && ~isempty(results_file)
                error('Cannot provide both results and results_file');
            elseif isempty(results) && isempty(results_file)
                error('Provide either results or results_file');
            elseif ~isempty(results)
                obj.results_cell = results;
                obj.results_table = obj.create_combined_table(results);
            else
                obj.load_results_from_file(results_file);
            end
        end
        
        function combined_table = create_combined_table(obj, results_cell)
            all_data = [];
            
            for i = 1:length(results_cell)
                result = results_cell{i};
                config = result.config;
                
                n_sims = length(result.num_arrivals);
                
                data_struct = struct();
                data_struct.C = repmat(config.C, n_sims, 1);
                data_struct.K = repmat(config.K, n_sims, 1);
                data_struct.lambda_arrival = repmat(config.lambda_arrival, n_sims, 1);
                data_struct.rho = repmat(config.rho, n_sims, 1);
                data_struct.num_simulations = repmat(config.num_simulations, n_sims, 1);
                data_struct.length_simulation = repmat(config.length_simulation, n_sims, 1);
                data_struct.num_arrivals = result.num_arrivals;
                data_struct.num_losses = result.num_losses;
                data_struct.average_num_users = result.average_num_users;
                data_struct.loss_probabilities = result.loss_probabilities;
                
                if isempty(all_data)
                    all_data = data_struct;
                else
                    fields = fieldnames(data_struct);
                    for j = 1:length(fields)
                        field = fields{j};
                        all_data.(field) = [all_data.(field); data_struct.(field)];
                    end
                end
            end
            
            combined_table = struct2table(all_data);
        end
        
        
        function grouped_stats = calculate_mean_and_group_by_ck(obj, column_name)
            if ~ismember(column_name, obj.results_table.Properties.VariableNames)
                error('Column %s not found in results table', column_name);
            end
            
            unique_combinations = unique(obj.results_table(:, {'C', 'K'}), 'rows');
            
            C_vals = unique_combinations.C;
            K_vals = unique_combinations.K;
            means = zeros(size(C_vals));
            stds = zeros(size(C_vals));
            counts = zeros(size(C_vals));
            
            for i = 1:length(C_vals)
                mask = obj.results_table.C == C_vals(i) & obj.results_table.K == K_vals(i);
                data = obj.results_table.(column_name)(mask);
                means(i) = mean(data);
                stds(i) = std(data);
                counts(i) = sum(mask);
            end
            
            grouped_stats = table(C_vals, K_vals, means, stds, counts, ...
                'VariableNames', {'C', 'K', 'Mean', 'Std', 'Count'});
            grouped_stats = sortrows(grouped_stats, {'C', 'K'});
        end
        
        function plot_mean_grouped_by_ck(obj, column_name, plot_title)
            if nargin < 3
                plot_title = strrep(column_name, '_', ' ');
                plot_title = title_case(plot_title);
            end
            
            grouped_stats = obj.calculate_mean_and_group_by_ck(column_name);
            
            figure('Position', [100, 100, 800, 600]);
            hold on;
            
            unique_C = unique(grouped_stats.C);
            colors = lines(length(unique_C));
            
            for i = 1:length(unique_C)
                c_val = unique_C(i);
                mask = grouped_stats.C == c_val;
                
                K_subset = grouped_stats.K(mask);
                means_subset = grouped_stats.Mean(mask);
                
                plot(K_subset, means_subset, '-o', 'Color', colors(i,:), ...
                    'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:), ...
                    'DisplayName', sprintf('C=%d', c_val));
            end
            
            xlabel('Buffer Length (K)', 'FontSize', 12);
            ylabel(plot_title, 'FontSize', 12);
            title(sprintf('%s by C and K', plot_title), 'FontSize', 14);
            grid on;
            legend('Location', 'best', 'FontSize', 11);
            hold off;
        end
    end
end