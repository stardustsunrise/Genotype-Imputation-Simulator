from __future__ import annotations
import numpy as np

MISSING = -1

def _reference(reference: np.ndarray) -> np.ndarray:
    return np.asarray(reference, dtype=np.int8)

def majority_impute(masked: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Impute every missing genotype using the modal reference genotype."""
    x = masked.copy()
    ref = _reference(reference)
    modes = np.array([
        np.bincount(ref[:, j], minlength=3).argmax()
        for j in range(ref.shape[1])
    ], dtype=np.int8)
    rows, cols = np.where(x == MISSING)
    x[rows, cols] = modes[cols]
    return x

def allele_frequency_impute(masked: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Draw an imputed genotype from Hardy-Weinberg genotype probabilities."""
    x = masked.copy()
    ref = _reference(reference)
    p = np.clip(ref.mean(axis=0) / 2.0, 1e-6, 1 - 1e-6)
    probs = np.column_stack(((1-p)**2, 2*p*(1-p), p**2))
    rng = np.random.default_rng(12345)
    rows, cols = np.where(x == MISSING)
    draws = np.array([
        rng.choice(3, p=probs[j])
        for j in cols
    ], dtype=np.int8)
    x[rows, cols] = draws
    return x

def ld_neighbor_impute(masked: np.ndarray, reference: np.ndarray,
                       window: int = 3) -> np.ndarray:
    """Simple LD-inspired imputer using nearby reference SNP correlation.

    This is deliberately a teaching approximation, not a production imputation
    algorithm. It predicts the missing dosage from the nearest correlated
    observed SNPs in the target individual.
    """
    x = masked.copy()
    ref = _reference(reference)

    for i in range(x.shape[0]):
        missing_cols = np.where(x[i] == MISSING)[0]
        for j in missing_cols:
            lo = max(0, j - window)
            hi = min(x.shape[1], j + window + 1)
            candidates = [k for k in range(lo, hi) if k != j and x[i, k] != MISSING]
            if not candidates:
                x[i, j] = majority_impute(x[i:i+1], ref)[i-i, j]
                continue

            best_k = None
            best_abs_r = -1.0
            target_k = candidates
            for k in target_k:
                a = ref[:, j].astype(float)
                b = ref[:, k].astype(float)
                if np.std(a) == 0 or np.std(b) == 0:
                    continue
                r = np.corrcoef(a, b)[0, 1]
                if np.isfinite(r) and abs(r) > best_abs_r:
                    best_abs_r = abs(r)
                    best_k = k

            if best_k is None:
                x[i, j] = majority_impute(x[i:i+1], ref)[0, j]
                continue

            # Simple nearest-neighbor dosage rule based on reference pairs.
            distances = np.abs(ref[:, best_k] - x[i, best_k])
            candidate_rows = np.argsort(distances)[:max(5, len(ref)//20)]
            predicted = int(np.rint(np.mean(ref[candidate_rows, j])))
            x[i, j] = np.clip(predicted, 0, 2)

    return x

def run_methods(masked: np.ndarray, reference: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "Majority genotype": majority_impute(masked, reference),
        "Allele-frequency": allele_frequency_impute(masked, reference),
        "LD-inspired": ld_neighbor_impute(masked, reference),
    }
