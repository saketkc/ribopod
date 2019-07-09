import pandas as pd
from plotly.graph_objs import Bar, Figure, Layout


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
        if column not in orf_counts_df.columns:
            continue
        if "translating" not in orf_counts_df[column].keys():
            continue
        trace = Bar(
            x=orf_counts_df.index,
            y=orf_counts_df[column]["translating"].tolist(),
            name=column,
        )
        fig_data.append(trace)
    layout = Layout(barmode="stack")
    fig = Figure(data=fig_data, layout=layout)
    return fig


def plot_phase_scores_violin(phase_score_df):
    """Create violoin plots for phase scores

    Parameters
    ----------
    phase_score_df: DataFrame
                    table withindex as ORF ID and columns as phase scores
    """

    columns = sorted(phase_score_df.columns.tolist())
    fig_data = []
    for column in columns:
        trace = {
            "type": "violin",
            "x": [column],
            "y": phase_score_df[column][phase_score_df[column] > 0],
            "name": column,
            "box": {"visible": True},
            "meanline": {"visible": True},
            # "line": {
            #    "color": 'black'
            # },
            # "points": 'all',
            # "jitter": 0,
        }
        fig_data.append(trace)
    layout = Layout({"title": "", "yaxis": {"zeroline": False}})
    fig = Figure(data=fig_data, layout=layout)
    return fig
