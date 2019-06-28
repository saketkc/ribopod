from collections import OrderedDict

import numpy as np
import pandas as pd
from plotly import tools
from plotly.graph_objs import Bar


def parse_ribotricer_bam_summary(file_path):
    """Parse ribotricer's bam_summary output.

    Parameters
    ----------
    file_path: string
               Location to bam_summar.txt

    Returns
    -------
    summary_dict: dict
                  Summary of mapped, unique reads
    fragment_length_dist: Series
                          A pandas series with index as fragment lengths
                          and values as read counts
    """

    summary_dict = OrderedDict()
    fragment_len_dist_dict = OrderedDict()
    reading_length_dist = False
    reading_summary = False
    with open(file_path) as fh:
        for index, line in enumerate(fh):
            line = line.strip()
            if line == "":
                continue
            if index == 0:
                assert line == "summary:"
                reading_summary = True
                continue
            if line == "length dist:":
                reading_summary = False
                reading_length_dist = True
                continue

            if reading_summary:
                try:
                    key, value = line.split(":")
                except:
                    raise Exception("Unable to parse {}".format(line))
                value = value.strip(" ")
                summary_dict[key] = int(value)
            if reading_length_dist:
                try:
                    key, value = line.split(":")
                except:
                    raise Exception("Unable to parse {}".format(line))
                value = value.strip(" ")
                fragment_len_dist_dict[int(key)] = int(value)
    return summary_dict, pd.Series(fragment_len_dist_dict).sort_index()


def project_summary_read_length_creator(project_summary_file):
    """Parse project summary file to prepare project for read length distribution analysis.

    Parameters
    ----------
    project_summary_file: string
                          Path to project summary file

    Returns
    -------
    read_length_dist_dict: dict
                           Keys as sample name, value as series returned by parse_ribotricer_bam_summary

    """
    summary_df = (
        pd.read_csv(project_summary_file, sep="\t")
        .set_index("experiment_accession")
        .sort_index()
    )
    summary_df = summary_df[
        summary_df.ribotricer_metagene_5p == summary_df.ribotricer_metagene_5p
    ]
    read_length_dist_dict = OrderedDict()
    for sample_name, row in summary_df.iterrows():
        # Load the 5' tsv
        if row.ribotricer_bam_summary:
            _, read_length_dist_dict[sample_name] = parse_ribotricer_bam_summary(
                row.ribotricer_bam_summary
            )
        else:
            read_length_dist_dict[sample_name] = pd.Series()
    return read_length_dist_dict


def plot_read_length_distribution(
    read_lengths, samples_per_row=1, shared_yaxes=False, plot_ridge=False
):
    """Plot read length distribution.

    Parameters
    ----------

    read_lengths : dict
                       Keys as sample name, value as series of read lengths
    """
    fig = tools.make_subplots(
        rows=int(np.ceil(len(list(read_lengths.keys())) / samples_per_row)),
        cols=samples_per_row,
        subplot_titles=list(read_lengths.keys()),
        print_grid=False,
    )

    index = 1
    if plot_ridge:
        pass
        # df = pd.DataFrame(read_lengths)
        # fig = draw_ridge_plot(df)

    else:
        for sample in list(read_lengths.keys()):
            trace = Bar(
                x=read_lengths[sample].index.tolist(),
                y=read_lengths[sample].sort_index().values.tolist(),
                name=sample,
            )
            fig.append_trace(
                trace,
                (index - 1) // samples_per_row + 1,
                (index - 1) % samples_per_row + 1,
            )
            index += 1
    fig["layout"].update(
        height=max(200 * len(list(read_lengths.keys())), 400),
        width=1000,
        title="Read length distribution",
    )
    # fig['layout'].update(scene=dict(aspectmode="data"))
    fig["layout"].update(font=dict(family="Arial", size=28, color="#000000"))
    fig["layout"].update(showlegend=False)
    return fig
