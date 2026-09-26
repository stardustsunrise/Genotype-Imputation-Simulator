import plotly.express as px
import pandas as pd

def metric_bar(results: pd.DataFrame, metric: str):
    return px.bar(
        results, x="Method", y=metric,
        title=f"{metric} by imputation method",
        range_y=[0, 1] if metric != "MAE" else None,
    )

def maf_plot(results: pd.DataFrame):
    return px.line(
        results, x="MAF bin", y="Accuracy", color="Method",
        markers=True, title="Accuracy across minor allele-frequency bins",
        range_y=[0, 1],
    )