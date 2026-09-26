import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import streamlit as st

from src.config import SimulationConfig
from src.simulation import simulate_genotypes
from src.imputation import run_methods
from src.evaluation import summarize_results, by_maf
from src.plots import metric_bar, maf_plot
from src.tool_catalog import TOOLS

st.set_page_config(page_title="Genotype Imputation Explorer", layout="wide")

st.title("🧬 Genotype Imputation Explorer")
st.caption("Educational simulator for understanding how imputation performance changes with data characteristics.")

with st.sidebar:
    st.header("Simulation")
    n_individuals = st.slider("Individuals", 50, 1000, 300, 50)
    n_snps = st.slider("SNPs", 20, 500, 100, 10)
    maf = st.slider("Minor allele frequency", 0.01, 0.50, 0.20, 0.01)
    ld = st.slider("LD strength", 0.0, 0.95, 0.70, 0.05)
    missing = st.slider("Missingness", 0.01, 0.50, 0.10, 0.01)
    ref_size = st.slider("Reference-panel size", 20, max(20, n_individuals - 1), min(150, n_individuals - 1), 10)
    seed = st.number_input("Random seed", 0, 100000, 42)

    run = st.button("Run simulation", type="primary", use_container_width=True)

if "result" not in st.session_state or run:
    config = SimulationConfig(
        n_individuals=n_individuals,
        n_snps=n_snps,
        maf=maf,
        ld_strength=ld,
        missing_rate=missing,
        reference_size=ref_size,
        seed=int(seed),
    )
    with st.spinner("Simulating genotypes and running methods..."):
        sim = simulate_genotypes(config)
        reference = sim.genotypes[sim.reference_indices]
        predictions = run_methods(sim.masked_genotypes, reference)
        results = summarize_results(sim.genotypes, predictions, sim.mask)
        maf_results = by_maf(sim.genotypes, predictions, sim.mask)
    st.session_state.result = (sim, predictions, results, maf_results)

sim, predictions, results, maf_results = st.session_state.result

st.subheader("Simulation summary")
c1, c2, c3, c4 = st.columns(4)
c1.metric("Individuals", sim.genotypes.shape[0])
c2.metric("SNPs", sim.genotypes.shape[1])
c3.metric("Masked genotypes", int(sim.mask.sum()))
c4.metric("Reference individuals", len(sim.reference_indices))

st.subheader("Performance")
st.dataframe(results.style.format({"Accuracy": "{:.3f}", "Dosage R²": "{:.3f}", "MAE": "{:.3f}"}), use_container_width=True)

left, right = st.columns(2)
with left:
    st.plotly_chart(metric_bar(results, "Accuracy"), use_container_width=True)
with right:
    st.plotly_chart(metric_bar(results, "Dosage R²"), use_container_width=True)

if not maf_results.empty:
    st.plotly_chart(maf_plot(maf_results), use_container_width=True)

st.subheader("Real-world imputation tools")
for tool in TOOLS:
    with st.expander(tool["name"]):
        st.write(f"**Category:** {tool['category']}")
        st.write(tool["description"])
        st.caption(tool["simulation_note"])

st.info(
    "The implemented methods are simplified educational baselines. "
    "They should not be treated as substitutes for production tools such as Beagle, "
    "IMPUTE, or Minimac."
)
