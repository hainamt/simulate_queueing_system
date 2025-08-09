classdef SimulationConfiguration
    %SIMULATIONCONFIGURATION Summary of this class goes here
    %   Detailed explanation goes here

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

        % function outputArg = method1(obj,inputArg)
        %     %METHOD1 Summary of this method goes here
        %     %   Detailed explanation goes here
        %     outputArg = obj.Property1 + inputArg;
        % end
    end
end