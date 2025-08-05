from typing import Optional

from configuration import SimulationConfiguration
from dataclasses import dataclass
import numpy as np

@dataclass
class VectorizedSimulationResult:
    config: SimulationConfiguration
    num_arrivals: np.ndarray
    num_losses: np.ndarray
    average_num_users: np.ndarray  # shape: (num_simulations,)
    loss_probabilities: np.ndarray  # shape: (num_simulations,)

    def to_dict(self):
        return {
            "config": self.config.to_dict(),
            "num_arrivals": self.num_arrivals.tolist(),
            "num_losses": self.num_losses.tolist(),
            "average_num_users": self.average_num_users.tolist(),
            "loss_probabilities": self.loss_probabilities.tolist(),
        }


class VectorizedQueueSimulator:

    def __init__(self, config: SimulationConfiguration):
        self.config = config
        self.C = config.C
        self.K = config.K

        self.lambda_arrival = config.lambda_arrival
        self.mu_service = config.mu_service

        self.num_arrival_stages = config.num_arrival_stages
        self.num_service_stages = config.num_service_stages

        self.num_simulations = config.num_simulations
        self.length_simulation = config.length_simulation

        self.num_arrivals = np.zeros(self.num_simulations, dtype=int)
        self.num_losses = np.zeros(self.num_simulations, dtype=int)

        self.x_in = np.ones(self.num_simulations, dtype=int)
        self.num_event_types = self.num_arrival_stages + self.num_service_stages * self.C
        self.event_table = np.full((self.num_simulations, self.num_event_types), np.inf)
        self.number_users_when_event = np.zeros((self.num_simulations, self.length_simulation))

        self.clk = np.zeros(self.num_simulations)
        self.previous_queue_change_time = np.zeros(self.num_simulations)

        self.total_customers = np.zeros(self.num_simulations, dtype=int)
        self.server_states = np.zeros((self.num_simulations, self.C), dtype=int)
        self.state_residence_time = np.zeros((self.num_simulations, self.K + 1))  # For states 0 to K

        self.active_simulations = np.ones(self.num_simulations, dtype=bool)

        base_seed = config.simulation_id
        self.rng_seeds = np.arange(base_seed, base_seed + self.num_simulations)
        self.rng = np.random.RandomState(base_seed)

        self.result: Optional[VectorizedSimulationResult] = None

    def rand_exp(self, sim_indices, is_service=True):
        size = len(sim_indices)
        rate = self.mu_service if is_service else self.lambda_arrival
        num_stages = self.num_service_stages if is_service else self.num_arrival_stages

        uniform_samples = self.rng.random(size)
        exponential_samples = -np.log(1 - uniform_samples) / (num_stages * rate)

        return exponential_samples if size > 1 else exponential_samples[0]

    def update_state_residence_time(self, sim_indices):
        time_diff = self.clk[sim_indices] - self.previous_queue_change_time[sim_indices]
        current_states = self.total_customers[sim_indices]

        self.state_residence_time[sim_indices, current_states] += time_diff

    def handle_a1_events(self, sim_indices):
        self.x_in[sim_indices] = 2

        self.event_table[sim_indices, 0] = np.inf
        self.event_table[sim_indices, 1] = self.clk[sim_indices] + self.rand_exp(sim_indices)

    def handle_a2_events(self, sim_indices):
        self.num_arrivals[sim_indices] += 1
        can_enter = self.total_customers[sim_indices] < self.K

        # get a list of simulations that have a customer entering the queue
        entering_sims = sim_indices[can_enter]

        if len(entering_sims) > 0:
            self.update_state_residence_time(entering_sims)
            self.total_customers[entering_sims] += 1

            # get server states simulation that has customers entering the queue
            server_states_entering = self.server_states[entering_sims]

            # for those simulations, get a mask of state servers
            idle_mask = server_states_entering == 0

            # check if there is an idle server in each simulation
            has_idle_server = np.any(idle_mask, axis=1)

            # create a list of assignable simulations (-1, -1, ... entering_sims)
            server_assignments = np.full(len(entering_sims), -1, dtype=int)

            if np.any(has_idle_server):
                # get index of first available server in each simulation
                first_idle_indices = np.argmax(idle_mask, axis=1)
                server_assignments[has_idle_server] = first_idle_indices[has_idle_server]

                # Filter to only simulations that actually have idle servers
                valid_mask = server_assignments >= 0
                valid_entering_sims = entering_sims[valid_mask]
                valid_server_assignments = server_assignments[valid_mask]

                if len(valid_entering_sims) > 0:
                    # Update server states vectorially
                    self.server_states[valid_entering_sims, valid_server_assignments] = 1

                    s1_event_indices = self.num_arrival_stages + 2 * valid_server_assignments
                    service_times = self.rand_exp(valid_entering_sims, is_service=True)
                    self.event_table[valid_entering_sims, s1_event_indices] = self.clk[valid_entering_sims] + service_times

            self.previous_queue_change_time[entering_sims] = self.clk[entering_sims]

        blocked = self.total_customers[sim_indices] >= self.K
        blocked_sims = sim_indices[blocked]
        if len(blocked_sims) > 0:
            self.num_losses[blocked_sims] += 1

        self.x_in[sim_indices] = 1
        self.event_table[sim_indices, 0] = self.clk[sim_indices] + self.rand_exp(sim_indices, is_service=False)
        self.event_table[sim_indices, 1] = np.inf

    def handle_s1_events(self,sim_indices, server_ids):
        self.server_states[sim_indices, server_ids] = 2
        s1_event_indices = 2 + 2 * server_ids
        s2_event_indices = s1_event_indices + 1

        self.event_table[sim_indices, s1_event_indices] = np.inf
        service_times = self.rand_exp(sim_indices, is_service=True)
        self.event_table[sim_indices, s2_event_indices] = self.clk[sim_indices] + service_times

    def handle_s2_events(self, sim_indices, server_ids):
        self.update_state_residence_time(sim_indices)
        self.total_customers[sim_indices] -= 1

        s1_event_indices = 2 + 2 * server_ids
        s2_event_indices = s1_event_indices + 1

        becomes_idle = self.total_customers[sim_indices] < self.C
        continues_serving = ~becomes_idle

        idle_sims = sim_indices[becomes_idle]
        idle_servers = server_ids[becomes_idle]
        if len(idle_sims) > 0:
            self.server_states[idle_sims, idle_servers] = 0
            idle_s1_indices = 2 + 2 * idle_servers
            self.event_table[idle_sims, idle_s1_indices] = np.inf

        serving_sims = sim_indices[continues_serving]
        serving_servers = server_ids[continues_serving]
        if len(serving_sims) > 0:
            self.server_states[serving_sims, serving_servers] = 1
            serving_s1_indices = 2 + 2 * serving_servers
            service_times = self.rand_exp(serving_sims)
            self.event_table[serving_sims, serving_s1_indices] = self.clk[serving_sims] + service_times

        self.event_table[sim_indices, s2_event_indices] = np.inf
        self.previous_queue_change_time[sim_indices] = self.clk[sim_indices]

    def simulate(self):
        active_sims = np.arange(self.num_simulations)
        self.event_table[active_sims, 0] = self.clk[active_sims] + self.rand_exp(active_sims, is_service=False)

        iteration = 0
        max_iterations = 1000000

        while np.any(self.active_simulations) and iteration < max_iterations:
            active_sims = np.where(self.active_simulations)[0]
            if len(active_sims) == 0:
                break

            # Find the next event for each active simulation
            next_event_times = np.min(self.event_table[active_sims], axis=1)
            next_event_indices = np.argmin(self.event_table[active_sims], axis=1)

            # Update clocks
            self.clk[active_sims] = next_event_times

            # Check which simulations should stop
            finished = self.clk[active_sims] >= self.config.length_simulation
            finished_sims = active_sims[finished]
            self.active_simulations[finished_sims] = False

            # Continue with still active simulations
            still_active = active_sims[~finished]
            if len(still_active) == 0:
                break

            still_active_event_indices = next_event_indices[~finished]

            # Group simulations by event type for vectorized processing
            for event_type in range(self.num_event_types):
                event_mask = still_active_event_indices == event_type
                if not np.any(event_mask):
                    continue

                event_sims = still_active[event_mask]

                if event_type == 0:  # A1 events
                    self.handle_a1_events(event_sims)
                elif event_type == 1:  # A2 events
                    self.handle_a2_events(event_sims)
                else:
                    server_id = (event_type - 2) // 2
                    if event_type % 2 == 0:  # S1 events (even indices >= 2)
                        server_ids = np.full(len(event_sims), server_id)
                        self.handle_s1_events(event_sims, server_ids)
                    else:  # S2 events (odd indices >= 3)
                        server_ids = np.full(len(event_sims), server_id)
                        self.handle_s2_events(event_sims, server_ids)

            iteration += 1

        self.calculate_final_statistics()

    def calculate_final_statistics(self):
        valid_arrivals = self.num_arrivals > 0
        loss_probabilities = np.zeros(self.num_simulations)
        loss_probabilities[valid_arrivals] = self.num_losses[valid_arrivals] / self.num_arrivals[valid_arrivals]

        total_times = np.sum(self.state_residence_time, axis=1)
        average_num_users = np.zeros(self.num_simulations)

        valid_times = total_times > 0
        if np.any(valid_times):
            states = np.arange(self.K + 1)
            weighted_sums = np.sum(self.state_residence_time * states[np.newaxis, :], axis=1)
            average_num_users[valid_times] = weighted_sums[valid_times] / total_times[valid_times]

        self.result = VectorizedSimulationResult(
            config=self.config,
            num_arrivals=self.num_arrivals.copy(),
            num_losses=self.num_losses.copy(),
            average_num_users=average_num_users.copy(),
            loss_probabilities=loss_probabilities.copy())


# def run_vectorized_simulations(sim_config: SimulationConfiguration) -> VectorizedQueueSimulator:
#     print("Running vectorized simulations with config:")
#     pprint(sim_config)
#     q_simulator = VectorizedQueueSimulator(sim_config)
#     start = time()
#     q_simulator.simulate()
#     end = time()
#     print(f"Simulation took {end - start} seconds")
#     return q_simulator
