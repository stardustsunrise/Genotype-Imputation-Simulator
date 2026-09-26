from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from .config import SimulationConfig

@dataclass
class SimulationResult:
    genotypes: np.ndarray
    masked_genotypes: np.ndarray
    mask: np.ndarray
    reference_indices: np.ndarray
    target_indices: np.ndarray

def _simulate_haplotype(n_haplotypes: int, n_snps: int, maf: float,
                        ld_strength: float, rng: np.random.Generator) -> np.ndarray:
    # Educational LD model: each SNP partly copies the previous SNP and
    # otherwise samples a new allele. This creates tunable local correlation.
    h = np.empty((n_haplotypes, n_snps), dtype=np.int8)
    h[:, 0] = rng.binomial(1, maf, size=n_haplotypes)
    for j in range(1, n_snps):
        copy = rng.random(n_haplotypes) < ld_strength
        fresh = rng.binomial(1, maf, size=n_haplotypes)
        h[:, j] = np.where(copy, h[:, j - 1], fresh)
    return h

def simulate_genotypes(config: SimulationConfig) -> SimulationResult:
    config.validate()
    rng = np.random.default_rng(config.seed)
    haplotypes = _simulate_haplotype(
        2 * config.n_individuals,
        config.n_snps,
        config.maf,
        config.ld_strength,
        rng,
    )
    genotypes = (haplotypes[0::2] + haplotypes[1::2]).astype(np.int8)

    indices = rng.permutation(config.n_individuals)
    reference_indices = np.sort(indices[:config.reference_size])
    target_indices = np.sort(indices[config.reference_size:])

    mask = np.zeros_like(genotypes, dtype=bool)
    mask[target_indices] = rng.random(
        (len(target_indices), config.n_snps)
    ) < config.missing_rate

    masked = genotypes.copy()
    masked[mask] = -1

    return SimulationResult(
        genotypes=genotypes,
        masked_genotypes=masked,
        mask=mask,
        reference_indices=reference_indices,
        target_indices=target_indices,
    )
