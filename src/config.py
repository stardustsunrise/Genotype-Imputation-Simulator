from dataclasses import dataclass

@dataclass
class SimulationConfig:
    n_individuals: int = 500
    n_snps: int = 100
    maf: float = 0.20
    ld_strength: float = 0.70
    missing_rate: float = 0.10
    reference_size: int = 250
    seed: int = 42

    def validate(self) -> None:
        if self.n_individuals < 10:
            raise ValueError("n_individuals must be at least 10")
        if self.n_snps < 2:
            raise ValueError("n_snps must be at least 2")
        if not 0 < self.maf <= 0.5:
            raise ValueError("maf must be in (0, 0.5]")
        if not 0 <= self.ld_strength <= 0.99:
            raise ValueError("ld_strength must be in [0, 0.99]")
        if not 0 < self.missing_rate < 1:
            raise ValueError("missing_rate must be in (0, 1)")
        if not 2 <= self.reference_size <= self.n_individuals:
            raise ValueError("reference_size must be between 2 and n_individuals")
