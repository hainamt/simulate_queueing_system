import numpy as np
import heapq
from typing import List, Optional
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class SimulationConfiguration:
    C: int
    K: int
    lambda_arrival: float
    mu_service: float
    num_arrival_stages: int
    num_service_stages: int
    length_simulation: int
    simulation_id: int


@dataclass
class SimulationResult:
    simulation_config: SimulationConfiguration
    num_arrivals: int
    num_losses: int
    average_num_users: float
    loss_probability: float


@dataclass
class Event:
    time: float
    event_type: str
    server_id: Optional[int] = None

    def __lt__(self, other):
        return self.time < other.time


class QueueSimulator:
    def __init__(self, config: SimulationConfiguration):
        self.config = config
        self.C = self.config.C
        self.K = self.config.K

        self.lambda_arrival = self.config.lambda_arrival
        self.mu_service = self.config.mu_service

        self.num_arrival_stages = self.config.num_arrival_stages
        self.num_service_stages = self.config.num_service_stages

        self.event_queue: List[Event] = []
        self.num_arrivals = 0
        self.num_losses = 0
        self.total_customers = 0
        self.number_users_when_event = dict()
        self.event_counter = 0

        self.state_residence_time = defaultdict(float)
        self.queue_change_time = 0

        self.servers_state = np.zeros(self.C)

        self.clk = 0

        self.acceptable_event_types = ['A' + str(i) for i in range(1, self.num_arrival_stages + 1)] + \
            ['S' + str(i) for i in range(1, self.num_service_stages + 1)]

    def update_state_residence_time(self):
        self.state_residence_time[self.total_customers] += (self.clk - self.queue_change_time)

    def schedule_event(self, event: Event):
        heapq.heappush(self.event_queue, event)

    def get_next_event(self) -> Event:
        return heapq.heappop(self.event_queue)

    def rand_exp_arrival(self):
        return np.random.exponential(1/(self.num_arrival_stages * self.lambda_arrival))

    def rand_exp_service(self):
        return np.random.exponential(1/(self.num_service_stages * self.mu_service))

    def handle_event(self, event: Event):

        timestamp = event.time
        event_type = event.event_type
        self.event_counter += 1
        self.clk = timestamp

        def handle_a1_event():
            next_event = Event(time=event.time + self.rand_exp_arrival(), event_type='A2')
            self.schedule_event(next_event)

        def handle_a2_event():
            self.num_arrivals += 1
            self.update_state_residence_time()
            if self.total_customers < self.K:
                self.total_customers += 1
                idle_servers = np.where(self.servers_state == 0)[0]
                if len(idle_servers) > 0:
                    server_id = int(idle_servers[0])
                    next_event = Event(time=event.time + self.rand_exp_service(),
                                       event_type='S1', server_id=server_id)
                    self.schedule_event(next_event)
                    self.servers_state[server_id] = 1

                self.queue_change_time = event.time
            else:
                self.num_losses += 1

            next_a1 = Event(time=event.time + self.rand_exp_arrival(), event_type='A1')
            self.schedule_event(next_a1)

        def handle_s1_event():
            server_id = event.server_id
            self.servers_state[server_id] = 2
            next_event = Event(time=event.time + self.rand_exp_service(),
                               event_type='S2', server_id=server_id)
            self.schedule_event(next_event)

        def handle_s2_event():
            server_id = event.server_id
            self.total_customers -= 1
            self.update_state_residence_time()
            if self.total_customers < self.C:
                self.servers_state[server_id] = 0
            else:
                self.servers_state[server_id] = 1
                next_event = Event(time=event.time + self.rand_exp_service(), event_type='S1', server_id=server_id)
                self.schedule_event(next_event)
            self.queue_change_time = event.time

        match event_type:
            case 'A1':
                handle_a1_event()
            case 'A2':
                handle_a2_event()
            case 'S1':
                handle_s1_event()
            case 'S2':
                handle_s2_event()
            case _:
                raise ValueError(f"Invalid event type: {event.event_type}, expected one of: {self.acceptable_event_types}")

        self.number_users_when_event[self.event_counter] = self.total_customers

    def simulate(self):
        first_event = Event(time=self.clk + self.rand_exp_arrival(), event_type='A1')
        self.schedule_event(first_event)
        while self.clk < self.config.length_simulation:
            next_event = self.get_next_event()
            self.handle_event(next_event)

        loss_probability = self.num_losses / self.num_arrivals if self.num_arrivals > 0 else 0
        print(f"Loss per arrival: {loss_probability}")

        total_time = sum(self.state_residence_time.values())
        if total_time > 0:
            average_num_users = sum(
                state * _time for state, _time in self.state_residence_time.items()
            ) / total_time
        else:
            average_num_users = 0

        return SimulationResult(
            simulation_config=self.config,
            num_arrivals=self.num_arrivals,
            num_losses=self.num_losses,
            average_num_users=average_num_users,
            loss_probability=loss_probability
        )


def run_single_simulation(config: SimulationConfiguration):
    # Set seed for reproducibility
    np.random.seed(config.simulation_id * 1000 + hash((config.C, config.K)) % 1000)

    simulator = QueueSimulator(config)
    result = simulator.simulate()
    return result


if __name__ == "__main__":
    pass
