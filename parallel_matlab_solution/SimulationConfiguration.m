classdef SimulationConfiguration
    properties
        simulation_id (1,1) uint16 = 1000
        lambda (1,1) double
        rho (1,1) double
        C (1,1) uint8 = 2
        K (1,1) uint8 = 2
        num_simulations (1,1) uint16 = 100
        length_simulation (1,1) double = 50000
        max_iterations (1,1) uint32 = 10
        mu (1,1) double
    end

    methods
        function obj = SimulationConfiguration(simulation_id, ...
                                                lambda, rho, ...
                                                C, K, ...
                                                num_sims, length_sim, max_iterations)
            obj.simulation_id = simulation_id;
            obj.lambda = lambda;
            obj.rho = rho;
            obj.C = C;
            obj.K = K;
            obj.num_simulations = num_sims;
            obj.length_simulation = length_sim;
            obj.mu = obj.lambda / (double(obj.C) * obj.rho);
            obj.max_iterations = max_iterations;
        end

        function dict = to_dict(obj)
            dict = struct();
            dict.simulation_id = obj.simulation_id;
            dict.lambda = obj.lambda;
            dict.mu = obj.mu;
            dict.rho = obj.rho;
            dict.C = obj.C;
            dict.K = obj.K;
            dict.num_simulations = obj.num_simulations;
            dict.length_simulation = obj.length_simulation;
        end
    end
end