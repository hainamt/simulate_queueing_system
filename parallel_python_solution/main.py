from configuration import SimulationConfiguration
from analysis import StatisticalAnalyzer, plot_confidence_intervals, plot_mean_grouped_by_ck
from simulation_utils import create_config_grid_ck, run_multiprocess_simulations
import os
import polars as pl
from time import time


if __name__ == "__main__":
    base_config = SimulationConfiguration(
        simulation_id=1000,
        lambda_arrival=1/10,
        rho=1.05,
        C=2,
        K=30,
        num_simulations=50,
        length_simulation=2000,
        max_iterations=1000000)

    # tic = time()
    # result =  run_single_simulation(base_config)
    # toc = time()
    # print(f"Time taken: {toc - tic:.2f} seconds")

    c_values = [2, 3]
    k_max = 30
    configs = create_config_grid_ck(base_config, c_values, k_max)
    print(f"Created {len(configs)} configurations")
    print(f"Maximum number of available processors: {os.cpu_count()}")
    start = time()
    results, results_file = run_multiprocess_simulations(
        configs=configs,
        num_processes=5,
        save_results=True,
        result_dir="./multiprocess_results")
    end = time()

    successful = sum(1 for result in results if result.error is None)
    failed = len(results) - successful

    print(f"Total time taken: {end - start:.2f} seconds")
    print(f"\nSimulation Complete!, Successful: {successful}, Failed: {failed}")

    confidence_level: float = 0.95
    analyzer = StatisticalAnalyzer(results)
    df_num_users_mean = analyzer.calculate_mean_and_group_by_ck(
        column=pl.col("average_num_users"),
        alias="mean_num_users")

    df_loss_prob_mean = analyzer.calculate_mean_and_group_by_ck(
        column=pl.col("loss_probabilities"),
        alias="mean_loss_prob")

    df_num_users_ci = analyzer.calculate_confidence_intervals_df(
        c_value=2,
        confidence_level=confidence_level)

    plot_confidence_intervals(df_num_users_ci)
    plot_mean_grouped_by_ck(df_num_users_mean, alias="mean_num_users")
    plot_mean_grouped_by_ck(df_loss_prob_mean,alias="mean_loss_prob")
