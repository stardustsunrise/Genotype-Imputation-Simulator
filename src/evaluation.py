from __future__ import annotations
import numpy as np
import pandas as pd

def genotype_accuracy(true: np.ndarray, predicted: np.ndarray, mask: np.ndarray) -> float:
    if not mask.any():
        return float("nan")
    return float(np.mean(true[mask] == predicted[mask]))

def dosage_r2(true: np.ndarray, predicted: np.ndarray, mask: np.ndarray) -> float:
    if mask.sum() < 2:
        return float("nan")
    y = true[mask].astype(float)
    yhat = predicted[mask].astype(float)
    if np.std(y) == 0 or np.std(yhat) == 0:
        return 0.0
    r = np.corrcoef(y, yhat)[0, 1]
    return float(r * r) if np.isfinite(r) else 0.0

def mae(true: np.ndarray, predicted: np.ndarray, mask: np.ndarray) -> float:
    if not mask.any():
        return float("nan")
    return float(np.mean(np.abs(true[mask].astype(float) - predicted[mask].astype(float))))

def summarize_results(true, predictions, mask) -> pd.DataFrame:
    rows = []
    for name, pred in predictions.items():
        rows.append({
            "Method": name,
            "Accuracy": genotype_accuracy(true, pred, mask),
            "Dosage R²": dosage_r2(true, pred, mask),
            "MAE": mae(true, pred, mask),
        })
    return pd.DataFrame(rows)

def by_maf(true: np.ndarray, predictions: dict, mask: np.ndarray,
           bins=(0.0, 0.05, 0.10, 0.20, 0.50)) -> pd.DataFrame:
    rows = []
    for j in range(true.shape[1]):
        allele_freq = true[:, j].mean() / 2
        maf = min(allele_freq, 1 - allele_freq)
        for name, pred in predictions.items():
            m = mask[:, j]
            if not m.any():
                continue
            rows.append({
                "SNP": j,
                "MAF": maf,
                "Method": name,
                "Accuracy": np.mean(true[m, j] == pred[m, j]),
            })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    labels = [f"{bins[i]:.2f}-{bins[i+1]:.2f}" for i in range(len(bins)-1)]
    out["MAF bin"] = pd.cut(out["MAF"], bins=bins, labels=labels, include_lowest=True)
    return out.groupby(["MAF bin", "Method"], observed=True)["Accuracy"].mean().reset_index()
