classdef E2E2CKSimulationResult
    properties
        config SimulationConfiguration
        num_arrivals (:,1) uint32
        num_losses (:,1) uint32
        average_num_users (:,1) double
        loss_probabilities (:,1) double
        state_residence_times (:,:) double
        state_probabilities (:,:) double
        total_simulation_times (:,1) double
        total_iterations (:,1) uint32
    end

    methods
        function obj = E2E2CKSimulationResult(config, num_arrivals, num_losses, ...
                average_num_users, loss_probabilities, state_residence_times, ...
                state_probabilities, total_simulation_times, total_iterations)
            
            obj.config = config;
            obj.num_arrivals = num_arrivals;
            obj.num_losses = num_losses;
            obj.average_num_users = average_num_users;
            obj.loss_probabilities = loss_probabilities;
            obj.state_residence_times = state_residence_times;
            obj.state_probabilities = state_probabilities;
            obj.total_simulation_times = total_simulation_times;
            obj.total_iterations = total_iterations; 
        end
        
        function dict = to_dict(obj)
            dict = struct();
            dict.config = obj.config.to_dict();
            dict.num_arrivals = obj.num_arrivals;
            dict.num_losses = obj.num_losses;
            dict.average_num_users = obj.average_num_users;
            dict.loss_probabilities = obj.loss_probabilities;
            dict.state_residence_times = obj.state_residence_times;
            dict.state_probabilities = obj.state_probabilities;
            dict.total_simulation_times = obj.total_simulation_times;
            dict.total_iterations = obj.total_iterations;
        end

        function table_result = to_table(obj)
            table_result = table(obj.num_arrivals, obj.num_losses, ...
                obj.average_num_users, obj.loss_probabilities, ...
                obj.total_simulation_times, ...
                'VariableNames', {'NumArrivals', 'NumLosses', 'AverageNumUsers', ...
                'LossProbabilities', 'TotalSimulationTimes'});
        end

    end
end