from collections import OrderedDict

import numpy as np
import pandas as pd
from plotly import tools
from plotly.graph_objs import Bar, Figure, Layout

import plotly.figure_factory as ff


def plot_orf_counts_stacked_bar(df):
    """Create stacked bar charts showing number of actively translating uORFs/ORFs/dORFs.

    Parameters
    ----------
    orf_counts_df: DataFrame
                table with index as the sample name and columns as # of uORF/dORF/
    """
    orf_counts_df = pd.pivot_table(
        df,
        columns=["ORF_type", "status"],
        index="experiment_accession",
        values=["count"],
    )[
        "count"
    ]  # .reset_index()

    columns = [
        "annotated",
        "super_uORF",
        "uORF",
        "overlap_uORF",
        "super_dORF",
        "dORF",
        "overlap_dORF",
        "novel",
    ]

    fig_data = []
    for column in columns:
        trace = Bar(
            x=orf_counts_df.index,
            y=orf_counts_df[column]["translating"].tolist(),
            name=column,
        )
        fig_data.append(trace)
    layout = Layout(barmode="stack")
    fig = Figure(data=fig_data, layout=layout)
    return fig
