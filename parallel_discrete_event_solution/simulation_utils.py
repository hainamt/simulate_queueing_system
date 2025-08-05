from vectorized_solution import VectorizedQueueSimulator
from configuration import SimulationConfiguration
from typing import List, Tuple
from datetime import datetime
import multiprocessing as mp
import os
import json


def create_config_grid_ck(base_config: SimulationConfiguration,
                          c_values: List[int], k_values: List[int]) -> List[SimulationConfiguration]:
    configs = []
    simulation_id = base_config.simulation_id
    for c in c_values:
        for k in k_values:
            if k >= c:
                config = SimulationConfiguration(
                    simulation_id=simulation_id,
                    lambda_arrival=base_config.lambda_arrival,
                    mu_service=base_config.mu_service,
                    num_arrival_stages=base_config.num_arrival_stages,
                    num_service_stages=base_config.num_service_stages,
                    C=c,
                    K=k,
                    num_simulations=base_config.num_simulations,
                    length_simulation=base_config.length_simulation)
                configs.append(config)
                simulation_id += 1
            else:
                print(f"Skipping invalid configuration: C={c}, K={k} (K must be >= C)")
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
                if "error" not in result_dict:
                    successful_count += 1
                    combined_result = {
                        "results": result_dict
                    }
                else:
                    failed_count += 1
                    combined_result = {
                        "error": result_dict["error"]
                    }

                f.write(json.dumps(combined_result, ensure_ascii=False) + "\n")

        print(f"All results saved to {results_file}")
        print(f"Total: {len(configs)}, Successful: {successful_count}, Failed: {failed_count}")

    return results
