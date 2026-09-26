import numpy as np
from src.config import SimulationConfig
from src.simulation import simulate_genotypes

def test_simulation_shapes():
    cfg = SimulationConfig(n_individuals=50, n_snps=20, reference_size=20)
    result = simulate_genotypes(cfg)
    assert result.genotypes.shape == (50, 20)
    assert result.masked_genotypes.shape == (50, 20)
    assert result.mask.shape == (50, 20)

def test_masked_values_are_missing_marker():
    cfg = SimulationConfig(n_individuals=50, n_snps=20, reference_size=20, missing_rate=0.2)
    result = simulate_genotypes(cfg)
    assert np.all(result.masked_genotypes[result.mask] == -1)
