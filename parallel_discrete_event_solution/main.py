from configuration import SimulationConfiguration
from parallel_discrete_event_solution.analysis import StatisticalAnalyzer
from simulation_utils import create_config_grid_ck, run_multiprocess_simulations
import os

if __name__ == "__main__":
    base_config = SimulationConfiguration(
        simulation_id=1000,
        lambda_arrival=1 / 10,
        mu_service=1 / (2 * 1.05),
        num_arrival_stages=2,
        num_service_stages=2,
        C=2,
        K=5,
        num_simulations=5,
        length_simulation=10000)

    c_values = [2, 3]
    k_max = 30

    configs = create_config_grid_ck(base_config, c_values, k_max)
    print(f"Created {len(configs)} configurations")
    print(f"Maximum number of processors: {os.cpu_count()}")
    results, results_file = run_multiprocess_simulations(
        configs=configs,
        num_processes=5,
        save_results=True,
        result_dir="./multiprocess_results")

    successful = sum(1 for result in results if "error" not in result)
    failed = len(results) - successful

    print(f"\nSimulation Complete!, Successful: {successful}, Failed: {failed}")

    # print("\nSample Results:")
    # for result in results[:3]:
    #     if "error" not in result:
    #         avg_loss = sum(result["loss_probabilities"]) / len(result["loss_probabilities"])
    #         avg_users = sum(result["average_num_users"]) / len(result["average_num_users"])
    #         print(f"  C={result['config']['C']}, K={result['config']['K']}: Loss={avg_loss:.4f}, Avg Users={avg_users:.2f}")

    confidence_level: float = 0.95
    analyzer = StatisticalAnalyzer(confidence_level)
    multiconfig_analysis = analyzer.analyze_multiple_results(results)
    analyzer.plot_results_with_ci()
