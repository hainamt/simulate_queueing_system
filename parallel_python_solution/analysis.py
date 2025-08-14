import polars as pl
import matplotlib.pyplot as plt
from typing import Dict, List, Optional
from parallel_python_solution.simulator import VectorizedE2E2CKSimulationResult

def plot_mean_grouped_by_ck(df: pl.DataFrame, alias: str):
    plt.figure(figsize=(10, 6))
    for c in df["C"].unique():
        mask = df["C"] == c
        plt.plot(df.filter(mask)["K"],
                 df.filter(mask)[alias],
                 label=f"C={c}",
                 marker='o')

    plt.xlabel("Buffer Length (K)")
    plt.ylabel(alias.title())
    plt.title(alias.title())
    plt.grid(True)
    plt.legend()
    plt.show()


class StatisticalAnalyzer:
    def __init__(self, confidence_level: float,
                 results: Optional[List[VectorizedE2E2CKSimulationResult]] = None,
                 results_file: Optional[str] = None):
        self.confidence_level = confidence_level

        if all([results, results_file]):
            raise ValueError("Can't provide both results and results_file")
        elif not any([results, results_file]):
            raise ValueError("Provide either results or results_file")
        elif results:
            df_list = [item.to_polar_df() for item in results]
            self.results_df = pl.concat(df_list, how="vertical")
        else:
            self.results_df = pl.read_ndjson(results_file)

        self.results_df = self.results_df.unnest("config")

    def calculate_mean_and_group_by_ck(self, column: pl.col, alias: str):
        df = (self.results_df.group_by("C", "K").agg(
            (column.list.eval(pl.element().sum()) / pl.col("num_simulations")).alias(
                alias)
        ).sort(["C", "K"]))
        return df.explode(alias).explode(alias)
