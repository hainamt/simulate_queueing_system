from dataclasses import dataclass, asdict


@dataclass
class SimulationConfiguration:
    simulation_id: int
    lambda_arrival: float
    mu_service: float
    num_arrival_stages: int
    num_service_stages: int
    C: int
    K: int

    num_simulations: int = 10000
    length_simulation: float = 50000

    @property
    def rho(self):
        return self.lambda_arrival / (self.C * self.mu_service)

    def to_dict(self):
        return asdict(self)
