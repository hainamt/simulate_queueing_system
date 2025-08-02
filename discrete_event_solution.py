import numpy as np
import heapq
from typing import List, Optional, Tuple, Dict
from dataclasses import dataclass
import concurrent.futures
from multiprocessing import cpu_count
import threading
import time


@dataclass
class Event:
    time: float
    event_type: str
    server_id: Optional[int] = None

    def __lt__(self, other):
        return self.time < other.time


@dataclass
class SimulationConfiguration:
    lambda_arrival: float
    mu_service: float
    num_arrival_stages: int
    num_service_stages: int
    num_servers_c: int
    C: int
    K: int

    @property
    def rho(self):
        return self.lambda_arrival / (self.C * self.mu_service)


class QueueSimulator:
    def __init__(self, configuration: SimulationConfiguration):
        self.event_queue : List[Event] = []
        self.C = configuration.C
        self.K = configuration.K

        self.lambda_arrival = configuration.lambda_arrival
        self.mu_service = configuration.mu_service

        self.num_arrival_stages = configuration.num_arrival_stages
        self.num_service_stages = configuration.num_service_stages

        self.event_queue: List[Event] = []
        self.num_arrivals = 0
        self.num_losses = 0
        self.total_customers = 0
        self.event_counter = 0

        self.servers_state = np.zeros(self.C)

        self.clk = 0

        self.acceptable_event_types = ['A' + str(i) for i in range(1, self.num_arrival_stages + 1)] + \
            ['S' + str(i) for i in range(1, self.num_service_stages + 1)]

    def schedule_event(self, event: Event):
        heapq.heappush(self.event_queue, event)

    def get_next_event(self) -> Event:
        return heapq.heappop(self.event_queue)

    def rand_exp_arrival(self):
        return np.random.exponential(1/(self.num_arrival_stages * self.lambda_arrival))

    def rand_exp_service(self):
        return np.random.exponential(1/(self.num_service_stages * self.mu_service))

    def handle_event(self, event: Event):
        def handle_a1_event():
            next_event = Event(time=event.time + self.rand_exp_arrival(), event_type='A2')
            self.schedule_event(next_event)

        def handle_a2_event():
            self.num_arrivals += 1
            if self.total_customers < self.K:
                self.total_customers += 1
                mask_server_state = self.servers_state == 0
                if np.any(mask_server_state):
                    server_id = int(np.argmin(mask_server_state))
                    next_event = Event(time=event.time + self.rand_exp_service(), event_type='S1', server_id=server_id)
                    self.schedule_event(next_event)
                    self.servers_state[server_id] = 1
            else:
                self.num_losses += 1

        def handle_s1_event():
            server_id = event.server_id
            self.servers_state[server_id] = 2
            next_event = Event(time=event.time + self.rand_exp_service(), event_type='S2', server_id=server_id)
            self.schedule_event(next_event)

        def handle_s2_event():
            server_id = event.server_id
            self.total_customers -= 1
            if self.total_customers <= self.C - 1:
                self.servers_state[server_id] = 0
            else:
                self.servers_state[server_id] = 1
                next_event = Event(time=event.time + self.rand_exp_service(), event_type='S1', server_id=server_id)
                self.schedule_event(next_event)

        timestamp = event.time
        event_type = event.event_type
        self.event_counter += 1
        self.clk = timestamp

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

    def simulate(self, simulation_duration: float = 50000):
        first_event = Event(time=self.clk + self.rand_exp_arrival(), event_type='A1')
        self.schedule_event(first_event)

        while self.clk < simulation_duration:
            next_event = self.get_next_event()
            self.event_counter += 1

            self.handle_event(next_event)



if __name__ == "__main__":
    configuration = SimulationConfiguration(
        lambda_arrival=10,
        mu_service=10,
        num_arrival_stages=2,
        num_service_stages=2,
        num_servers_c=2,
        C=2,
        K=2
    )







