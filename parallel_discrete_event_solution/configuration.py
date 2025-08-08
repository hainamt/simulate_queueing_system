from dataclasses import dataclass, asdict


@dataclass
class SimulationConfiguration:
    simulation_id: int
    lambda_arrival: float
    rho: float
    num_arrival_stages: int
    num_service_stages: int
    C: int
    K: int

    num_simulations: int = 10000
    length_simulation: float = 50000

    @property
    def mu_service(self):
        return self.lambda_arrival / (self.C * self.rho)

    def to_dict(self):
        return asdict(self)

