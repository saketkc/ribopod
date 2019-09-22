import pandas as pd
from plotly.graph_objs import Bar, Figure, Layout


def format_figure(fig):
    fig.update_xaxes(
        showgrid=False,
        showline=True,
        linewidth=2,
        linecolor="black",
        gridwidth=1,
        gridcolor="Gray",
        tickwidth=1,
        ticklen=10,
        ticks="outside",
        tickcolor="black",
    )
    fig.update_yaxes(
        showgrid=True,
        showline=True,
        linewidth=2,
        linecolor="black",
        gridwidth=1,
        gridcolor="Gray",
        tickwidth=1,
        ticklen=10,
        ticks="outside",
        tickcolor="black",
    )


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
    layout = Layout(
        barmode="stack", paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)"
    )
    fig = Figure(data=fig_data, layout=layout)
    format_figure(fig)
    fig["layout"].update(font=dict(family="Arial", size=18, color="#000000"))
    return fig


def plot_phase_scores_violin(phase_score_df):
    """Create violoin plots for phase scores

    Parameters
    ----------
    phase_score_df: DataFrame
                    table withindex as ORF ID and columns as phase scores
    """
    if phase_score_df is None:
        return None
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
    layout = Layout(
        {
            "title": "",
            "yaxis": {"zeroline": False},
            "paper_bgcolor": "rgba(0,0,0,0)",
            "plot_bgcolor": "rgba(0,0,0,0)",
        }
    )
    fig = Figure(data=fig_data, layout=layout)
    format_figure(fig)
    fig["layout"].update(font=dict(family="Arial", size=18, color="#000000"))
    return fig
