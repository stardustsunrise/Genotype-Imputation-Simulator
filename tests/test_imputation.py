import numpy as np
from src.imputation import run_methods

def test_methods_fill_missing_values():
    rng = np.random.default_rng(1)
    true = rng.integers(0, 3, size=(30, 10), dtype=np.int8)
    masked = true.copy()
    mask = rng.random(true.shape) < 0.2
    masked[mask] = -1
    reference = true[:15]
    outputs = run_methods(masked, reference)
    for pred in outputs.values():
        assert np.all(pred >= 0)
        assert np.all(pred <= 2)
        assert pred.shape == true.shape
