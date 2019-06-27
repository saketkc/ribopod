from collections import OrderedDict

import numpy as np
import pandas as pd
from plotly import tools
from plotly.graph_objs import Scatter, Bar  # , Heatmap Figure,

# import plotly.figure_factory as ff

__FRAME_COLORS__ = ["#fc8d62", "#66c2a5", "#8da0cb"]


def parse_metagene_profile(file_path):
    """Parse metagene profile file

    Parameters
    ----------
    file_path: string
               Path to metagene 5p profile

    Returns
    -------
    fragment_length_profile: dict
                             An OrderedDict of

    """
    # fragment_length	offset_5p	profile	phase_score	valid_codons
    # Set the index on fragment length for easy retrieval
    metagene_df = pd.read_csv(file_path, sep="\t").set_index("fragment_length")
    metagene_df["profile"] = [eval(profile) for profile in metagene_df.profile]
    metagene_df["profile"] = [
        pd.Series(profile, index=range(-offset_5p, len(profile) - offset_5p))
        for offset_5p, profile in zip(metagene_df.offset_5p, metagene_df.profile)
    ]
    return metagene_df


def project_summary_metagene_creator(project_summary_file):
    """Parse project summary file to prepare project for metagene analysis.

    Parameters
    ----------
    project_summary_file: string
                          Path to project summary file

    Returns
    -------
    metagene_dfs: dict
                  Keys as sample name, value as df loaded through parse_metagene_profile_file

    """
    summary_df = (
        pd.read_csv(project_summary_file, sep="\t")
        .set_index("experiment_accession")
        .sort_index()
    )
    summary_df = summary_df[summary_df.ribotricer_metagene_5p==summary_df.ribotricer_metagene_5p]
    metagene_dfs = OrderedDict()
    for sample_name, row in summary_df.iterrows():
        # Load the 5' tsv
        if row.ribotricer_metagene_5p:
            metagene_dfs[sample_name] = parse_metagene_profile(
                row.ribotricer_metagene_5p
            )
        else:
            metagene_dfs[sample_name] = pd.DataFrame()
    return metagene_dfs


def plot_metagene_coverage(
    metagene_dfs,
    fragment_length,
    position_range="default",
    samples_per_row=1,
    plot_type="bar",
    normalize_per_codon=False,
    shared_yaxes=False,
):
    """Plot metagene for a list of samples.

    Parameters
    ----------

    metagene_dfs: dict
                  Keys as sample name, value as df loaded through parse_metagene_profile_file
    fragment_length: int
                     Fragment length for plotting
    position_range: range
                    Range of position to plot
    samples_per_row: int
                     Number of samples to plot
    plot_type: string
               Options of bar or line
    """
    titles = list(metagene_dfs.keys())

    fig = tools.make_subplots(
        rows=int(np.ceil(len(list(metagene_dfs.keys())) / samples_per_row)),
        cols=samples_per_row,
        subplot_titles=titles,
        print_grid=False,
    )

    index = 1
    for sample_name, metagene_df in metagene_dfs.items():
        # Pull out the profile of coverage
        try:
            row = metagene_df.loc[fragment_length]
        except KeyError:
            continue

        profile = row.profile
        if position_range != "default":
            profile = profile[position_range]
        else:
            position_range = profile.index
        phase_score = row.phase_score

        if plot_type == "bar":
            frame0 = profile[profile.index % 3 == 0]
            frame1 = profile[profile.index % 3 == 1]
            frame2 = profile[profile.index % 3 == 2]
            trace0 = Bar(
                x=frame0.index,
                y=frame0.values,
                name="Frame 0",
                marker=dict(color=__FRAME_COLORS__[0]),
            )
            trace1 = Bar(
                x=frame1.index,
                y=frame1.values,
                name="Frame 1",
                marker=dict(color=__FRAME_COLORS__[1]),
            )
            trace2 = Bar(
                x=frame2.index,
                y=frame2.values,
                name="Frame 2",
                marker=dict(color=__FRAME_COLORS__[2]),
            )
            trace = [trace0, trace1, trace2]
            fig.append_trace(
                trace0,
                (index - 1) // samples_per_row + 1,
                (index - 1) % samples_per_row + 1,
            )

            fig.append_trace(
                trace1,
                (index - 1) // samples_per_row + 1,
                (index - 1) % samples_per_row + 1,
            )
            fig.append_trace(
                trace2,
                (index - 1) // samples_per_row + 1,
                (index - 1) % samples_per_row + 1,
            )
        elif plot_type == "line":
            trace = Scatter(x=position_range, y=profile.values, mode="lines")
            fig.append_trace(
                trace,
                (index - 1) // samples_per_row + 1,
                (index - 1) % samples_per_row + 1,
            )
        else:
            raise Exception(
                "Unknown plot_type '{}' selected. Possible options: 'bar' or 'line'.".format(
                    plot_type
                )
            )
        index += 1
    fig["layout"].update(
        height=max(200 * len(list(metagene_dfs.keys())), 400),
        width=1000,
        title="Metagene distribution",
    )
    fig["layout"].update(font=dict(family="Arial", size=28, color="#000000"))
    for i in fig["layout"]["annotations"]:
        i["font"] = dict(size=28)  # ,color='#ff0000')
    fig["layout"].update(showlegend=False)
    return fig
