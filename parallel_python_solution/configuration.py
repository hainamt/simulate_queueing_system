from dataclasses import dataclass, asdict


@dataclass
class SimulationConfiguration:
    simulation_id: int
    lambda_arrival: float
    rho: float
    C: int
    K: int

    num_simulations: int = 100
    length_simulation: float = 50000

    max_iterations: int = 1000000

    @property
    def mu_service(self):
        return self.lambda_arrival / (self.C * self.rho)

    def to_dict(self):
        return asdict(self)
