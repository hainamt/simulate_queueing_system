import polars as pl
import numpy as np
from scipy.stats import t as stats_t
import matplotlib.pyplot as plt
from typing import List, Optional
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


def plot_confidence_intervals(ci_df: pl.DataFrame, print_statistics: bool = True) -> None:
    k_values = ci_df["K"].to_numpy()
    means = ci_df["mean"].to_numpy()
    lower_ci = ci_df["ci_lower"].to_numpy()
    upper_ci = ci_df["ci_upper"].to_numpy()

    lower_errors = means - lower_ci
    upper_errors = upper_ci - means

    plt.figure(figsize=(12, 8))
    plt.errorbar(k_values, means,
                 yerr=[lower_errors, upper_errors],
                 fmt='o-',
                 capsize=5,
                 capthick=2,
                 linewidth=2,
                 markersize=6)

    plt.xlabel("K", fontsize=12)
    plt.ylabel("Average number of users", fontsize=12)
    plt.title(f"Average number of users (CI) - C = {ci_df['C'][0]}", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.xticks(k_values)
    y_min = np.floor(lower_ci.min())
    y_max = np.ceil(upper_ci.max())
    plt.yticks(np.arange(y_min, y_max + 1, 2))
    plt.tight_layout()
    plt.show()

    # if print_statistics:
    #     print(f"\nConfidence Interval Statistics for C = {ci_df['C'][0]}:")
    #     print(f"Confidence Level: {ci_df['confidence_level'][0] * 100}%")
    #     print("\nK\tMean\t\tCI Lower\tCI Upper\tWidth\t\tStd Dev\t\tn")
    #     print("-" * 80)
    #     for row in ci_df.iter_rows(named=True):
    #         width = row["ci_upper"] - row["ci_lower"]
    #         print(f"{row['K']}\t{row['mean']:.4f}\t\t{row['ci_lower']:.4f}\t\t"
    #               f"{row['ci_upper']:.4f}\t\t{width:.4f}\t\t"
    #               f"{row['std']:.4f}\t\t{row['n']}")


class StatisticalAnalyzer:
    def __init__(self,
                 results: Optional[List[VectorizedE2E2CKSimulationResult]] = None,
                 results_file: Optional[str] = None):

        if all([results, results_file]):
            raise ValueError("Can't provide both results and results_file")
        elif not any([results, results_file]):
            raise ValueError("Provide either results or results_file")
        elif results:
            df_list = [item.to_df_dict() for item in results]
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

    def calculate_confidence_intervals_df(
            self,
            c_value: int,
            confidence_level: float = 0.95):
        c_filtered_df = self.results_df.filter(pl.col("C") == c_value)
        if c_filtered_df.is_empty():
            raise ValueError(f"No data found for C = {c_value}")

        exploded_df = c_filtered_df.explode("average_num_users")

        grouped = (
            exploded_df
            .group_by("K")
            .agg([
                pl.col("average_num_users").mean().alias("mean"),
                pl.col("average_num_users").std(ddof=1).alias("std"),
                pl.count("average_num_users").alias("n")
            ])
            .sort("K")
        )

        alpha = (1 - confidence_level) / 2
        dfree = (grouped["n"] - 1).to_numpy()
        t_critical = stats_t.ppf(1 - alpha, dfree)
        margin_error = t_critical * (grouped["std"] / np.sqrt(grouped["n"]))

        return grouped.with_columns([
            (grouped["mean"] - margin_error).alias("ci_lower"),
            (grouped["mean"] + margin_error).alias("ci_upper"),
            pl.lit(confidence_level).alias("confidence_level"),
            pl.lit(c_value).alias("C")
        ])



