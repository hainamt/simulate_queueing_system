from vectorized_solution import VectorizedQueueSimulator
from configuration import SimulationConfiguration
from typing import List
from datetime import datetime
import multiprocessing as mp
import os
import json


def create_config_grid_ck(
        base_config: SimulationConfiguration,
        c_values=None,
        k_max: int = 30
) -> List[SimulationConfiguration]:
    if c_values is None:
        c_values = [2, 3]
    configs = []
    simulation_id = base_config.simulation_id

    for c in c_values:
        k_values = list(range(c, k_max + 1))
        for k in k_values:
            config = SimulationConfiguration(
                simulation_id=simulation_id,
                lambda_arrival=base_config.lambda_arrival,
                mu_service=base_config.mu_service,
                num_arrival_stages=base_config.num_arrival_stages,
                num_service_stages=base_config.num_service_stages,
                C=c,
                K=k,
                num_simulations=base_config.num_simulations,
                length_simulation=base_config.length_simulation
            )
            configs.append(config)
            simulation_id += 1

    return configs

def run_single_simulation(config: SimulationConfiguration):
    try:
        print(f"Starting simulation: C={config.C}, K={config.K}, ID={config.simulation_id}")
        simulator = VectorizedQueueSimulator(config)
        simulator.simulate()
        result_dict = simulator.result.to_dict()
        print(f"Completed simulation: C={config.C}, K={config.K}, ID={config.simulation_id}")
        return result_dict
    except Exception as e:
        print(f"Error in simulation C={config.C}, K={config.K}, ID={config.simulation_id}: {str(e)}")
        return {"error": str(e)}


def run_multiprocess_simulations(
        configs: List[SimulationConfiguration],
        num_processes: int = None,
        save_results: bool = True,
        result_dir: str = "./results"):
    if num_processes is None:
        num_processes = min(mp.cpu_count(), len(configs))

    print(f"Running {len(configs)} simulations using {num_processes} processes")
    if save_results:
        os.makedirs(result_dir, exist_ok=True)

    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(run_single_simulation, configs)

    if save_results:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        results_file = os.path.join(result_dir, f"{timestamp}_results.jsonl")
        successful_count = 0
        failed_count = 0

        with open(results_file, "w", encoding="utf-8") as f:
            for result_dict in results:
                f.write(json.dumps(result_dict, ensure_ascii=False) + "\n")

        print(f"All results saved to {results_file}")
        print(f"Total: {len(configs)}, Successful: {successful_count}, Failed: {failed_count}")

    return results
