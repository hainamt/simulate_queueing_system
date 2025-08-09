classdef SimulationConfiguration
    properties
        simulation_id (1,1) uint16 = 1000
        lambda (1,1) double {mustBePositive}
        mu (1,1) double {mustBePositive}
        rho (1,1) double {mustBePositive}
        num_arrival_stages (1,1) uint8 = 2
        num_service_stages (1,1) uint8 = 2
        C (1,1) uint8 = 2
        K (1,1) uint8 = 2
        num_simulations (1,1) uint16 = 100
        length_simulation (1,1) double {mustBePositive} = 50000
    end

    methods
        function obj = SimulationConfiguration(simulation_id, ...
                                                lambda, rho, ...
                                                num_arrival_stages, num_service_stages, ...
                                                C, K, ...
                                                num_sims, length_sim)
            obj.simulation_id = simulation_id;
            obj.lambda = lambda;
            obj.rho = rho;
            obj.num_arrival_stages = num_arrival_stages;
            obj.num_service_stages = num_service_stages;
            obj.C = C;
            obj.K = K;
            obj.num_simulations = num_sims;
            obj.length_simulation = length_sim;
            obj.mu = obj.lambda / (obj.C * obj.rho);
        end

        function dict = to_dict(obj)
            dict = struct();
            dict.simulation_id = obj.simulation_id;
            dict.lambda = obj.lambda;
            dict.mu = obj.mu;
            dict.rho = obj.rho;
            dict.num_arrival_stages = obj.num_arrival_stages;
            dict.num_service_stages = obj.num_service_stages;
            dict.C = obj.C;
            dict.K = obj.K;
            dict.num_simulations = obj.num_simulations;
            dict.length_simulation = obj.length_simulation;
        end
    end
end