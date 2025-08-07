import numpy as np
from scipy import stats
from dataclasses import dataclass
from typing import Dict, List, Optional
import matplotlib.pyplot as plt
from configuration import SimulationConfiguration


@dataclass
class StatisticalSummary:
    config: SimulationConfiguration

    loss_prob_mean: float
    loss_prob_var: float
    loss_prob_std: float
    loss_prob_ci_lower_upper: List[float]
    loss_prob_ci_width: float
    loss_prob_relative_width: float

    avg_users_mean: float
    avg_users_var: float
    avg_users_std: float
    avg_users_ci_lower_upper: List[float]
    avg_users_ci_width: float
    avg_users_relative_width: float

    state_probabilities: np.ndarray
    state_prob_std: np.ndarray

    confidence_level: float
    degrees_freedom: int
    t_critical: float


@dataclass
class MultiConfigurationAnalysis:
    summaries: List[StatisticalSummary]
    c_values: List[int]
    k_values: List[int]

    def get_summary_by_config(self, c: int, k: int) -> Optional[StatisticalSummary]:
        for summary in self.summaries:
            if summary.config.C == c and summary.config.K == k:
                return summary
        return None


class StatisticalAnalyzer:
    def __init__(self, confidence_level: float = 0.95):
        self.confidence_level = confidence_level
        self.alpha = (1 - confidence_level) / 2
        self.analysis: Optional[MultiConfigurationAnalysis] = None

    def analyze_single_result(self, result_dict: Dict) -> StatisticalSummary:
        config_dict = result_dict['config']
        config = SimulationConfiguration(**config_dict)

        loss_probs = np.array(result_dict['loss_probabilities'])
        avg_users = np.array(result_dict['average_num_users'])

        num_sims = len(loss_probs)
        df = num_sims - 1
        t_crit = stats.t.ppf(1 - self.alpha, df) if df > 0 else np.inf

        def calculate_statistics(arr):
            mean = np.mean(arr)
            var = np.var(arr, ddof=1)
            std = np.sqrt(var)
            sem = std / np.sqrt(num_sims)
            ci_half = t_crit * sem
            ci_lower = mean - ci_half
            ci_upper = mean + ci_half
            ci_width = ci_upper - ci_lower
            relative_width = ci_width / mean if mean > 0 else np.inf
            return mean, var, std, ci_lower, ci_upper, ci_width, relative_width

        (loss_mean, loss_var, loss_std,
         loss_ci_lower, loss_ci_upper, loss_ci_width, loss_relative_width) = calculate_statistics(loss_probs)

        (users_mean, users_var, users_std,
         users_ci_lower, users_ci_upper, users_ci_width, users_relative_width) = calculate_statistics(avg_users)

        state_probs = np.zeros(config.K + 1)
        state_prob_std = np.zeros(config.K + 1)

        return StatisticalSummary(
            config=config,
            loss_prob_mean=loss_mean,
            loss_prob_var=loss_var,
            loss_prob_std=loss_std,
            loss_prob_ci_lower_upper=[loss_ci_lower, loss_ci_upper],
            loss_prob_ci_width=loss_ci_width,
            loss_prob_relative_width=loss_relative_width,
            avg_users_mean=users_mean,
            avg_users_var=users_var,
            avg_users_std=users_std,
            avg_users_ci_lower_upper = [users_ci_lower, users_ci_upper],
            avg_users_ci_width=users_ci_width,
            avg_users_relative_width=users_relative_width,
            state_probabilities=state_probs,
            state_prob_std=state_prob_std,
            confidence_level=self.confidence_level,
            degrees_freedom=df,
            t_critical=t_crit
        )

    def analyze_multiple_results(self, results: List[Dict]):
        summaries = []
        c_values = set()
        k_values = set()

        for result_dict in results:
            if 'error' not in result_dict:
                summary = self.analyze_single_result(result_dict)
                summaries.append(summary)
                c_values.add(summary.config.C)
                k_values.add(summary.config.K)

        self.analysis = MultiConfigurationAnalysis(
            summaries=summaries,
            c_values=sorted(list(c_values)),
            k_values=sorted(list(k_values))
        )

    def create_comparison_table(self, analysis: MultiConfigurationAnalysis) -> str:
        lines = ["Configuration Analysis Summary", "=" * 80,
                 f"Confidence Level: {self.confidence_level * 100:.1f}%",
                 ""]
        header = f"{'C':>3} {'K':>3} {'Loss Mean':>10} {'Loss CI':>15} {'Users Mean':>11} {'Users CI':>15} {'N':>5}"
        lines.append(header)
        lines.append("-" * 80)

        for summary in analysis.summaries:
            c = summary.config.C
            k = summary.config.K
            loss_mean = summary.loss_prob_mean
            loss_ci = f"±{summary.loss_prob_ci_width / 2:.4f}"
            users_mean = summary.avg_users_mean
            users_ci = f"±{summary.avg_users_ci_width / 2:.2f}"
            n = summary.config.num_simulations

            row = f"{c:>3} {k:>3} {loss_mean:>10.4f} {loss_ci:>15} {users_mean:>11.2f} {users_ci:>15} {n:>5}"
            lines.append(row)

        return "\n".join(lines)


    def plot_results_with_ci(self, save_path: Optional[str] = None):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        for c in self.analysis.c_values:
            c_summaries = [s for s in self.analysis.summaries if s.config.C == c]
            c_summaries.sort(key=lambda x: x.config.K)

            k_vals = [s.config.K for s in c_summaries]

            users_means = [s.avg_users_mean for s in c_summaries]
            users_errors = [s.avg_users_ci_width / 2 for s in c_summaries]

            ax1.errorbar(k_vals, users_means, yerr=users_errors,
                         label=f'C={c}', linewidth=2, capsize=5, capthick=2)

            loss_means = [s.loss_prob_mean for s in c_summaries]
            loss_errors = [s.loss_prob_ci_width / 2 for s in c_summaries]

            ax2.errorbar(k_vals, loss_means, yerr=loss_errors,
                         label=f'C={c}', linewidth=2, capsize=5, capthick=2)

        ax1.set_title('Average Number of Users with Confidence Intervals')
        ax1.set_xlabel('K')
        ax1.set_ylabel('Average Number of Users')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_xticks(range(min(self.analysis.k_values), max(self.analysis.k_values) + 1))

        ax2.set_title('Loss Probability with Confidence Intervals')
        ax2.set_xlabel('K')
        ax2.set_ylabel('Loss Probability')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_xticks(range(min(self.analysis.k_values), max(self.analysis.k_values) + 1))

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        else:
            plt.show()

        return fig, (ax1, ax2)


def load_and_analyze_results(results_file: str, confidence_level: float = 0.95):
    import json
    results = []
    with open(results_file, 'r') as f:
        for line in f:
            if line.strip():
                results.append(json.loads(line))

    analyzer = StatisticalAnalyzer(confidence_level)
    return analyzer.analyze_multiple_results(results)
