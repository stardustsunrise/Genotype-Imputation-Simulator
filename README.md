# Genotype Imputation Explorer

A small educational/research-oriented Python project for simulating genotype data,
masking observed genotypes, applying simplified imputation strategies, and comparing
their performance.

This project intentionally does **not** reimplement production tools such as Beagle,
IMPUTE5, or Minimac. It provides simplified methods for studying the statistical ideas
behind genotype imputation.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

## Project structure

- `app.py` - Streamlit interface
- `src/simulation.py` - genotype and LD simulation
- `src/imputation.py` - simplified imputation methods
- `src/evaluation.py` - accuracy, dosage R², MAE, and stratified metrics
- `src/config.py` - dataclasses and defaults
- `src/tool_catalog.py` - metadata for real-world imputation tools
- `src/plots.py` - reusable Plotly figures
- `tests/` - basic tests
- `data/` - optional local datasets
- `notebooks/` - optional exploratory notebooks

## Scientific scope

The simulator uses diploid genotype dosages coded 0/1/2. It can generate correlated
variants using a simple latent haplotype process, mask genotypes, and compare several
educational baselines.

For a research-grade project, validate assumptions and compare against established
software using appropriate benchmark datasets.
