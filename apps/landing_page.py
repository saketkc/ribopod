import dash_core_components as dcc
import dash_html_components as html

import dash_table as dt
from dash.dependencies import Input, Output

from dash.dependencies import State
from init import __DATASETS__
from init import __SPECIES__

from metagene_helper import project_summary_metagene_creator, plot_metagene_coverage
from project_helper import (
    get_projects,
    get_srp_table,
    get_srp_read_lengths,
    get_project_summary_file,
)

from dash_helper import generate_table
from app_multipage import app


layout = html.Div(
    [
        html.Div(
            [
                html.H2(
                    "ribopod: A database of actively-translating ORFs in Ribo-seq data"
                ),
                html.Hr(),
            ]
        ),
        html.Div(id="srametadata-table"),
        html.Div(id="datatable-interactivity-container"),
        html.Div(
            [
                html.Div(
                    [
                        html.Label("Species: "),
                        dcc.Dropdown(options=__SPECIES__, value="hg38", id="assembly"),
                    ],
                    style={
                        "width": "25%",
                        "display": "inline-block",
                        "font-family": "Droid Serif",
                    },
                ),
                html.Div(
                    [
                        html.Label("SRP: "),
                        dcc.Dropdown(options={}, value=None, id="srp"),
                    ],
                    style={
                        "width": "25%",
                        "display": "inline-block",
                        "font-family": "Droid Serif",
                    },
                ),
                html.Div(
                    [
                        html.Label("Gene: "),
                        dcc.Input(
                            id="gene-input", type="text", placeholder="Gene name"
                        ),
                        html.Button("Submit", id="gene-submit"),
                    ],
                    style={
                        "width": "25%",
                        "font-family": "Droid Serif",
                        "display": "none",
                    },
                ),
            ]
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.Label("Read length for metagene plot: "),
                                dcc.Dropdown(
                                    options=28, value=28, id="read_length"
                                ),
                            ],
                            style={
                                "width": "25%",
                                "display": "inline-block",
                                "font-family": "Droid Serif",
                            },
                        ),
                        html.Div(
                            [
                                dcc.Checklist(
                                    options=[
                                        {
                                            "label": "Codon level normalization (automatically de-trends)",
                                            "value": True,
                                        }
                                    ],
                                    id="normalize-chk",
                                    values=[False],
                                )
                            ],
                            style={
                                "width": "65%",
                                "display": "inline-block",
                                "height": "100%",
                                "align": "center",
                            },
                        ),
                    ]
                )
            ]
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [dcc.Graph(id="coherence-heatmap")],
                            style={
                                "width": "100%",
                                "display": "inline-block",
                                "height": "100%",
                                "align": "center",
                            },
                        )
                    ]
                )
            ]
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [dcc.Graph(id="length-dist-plot")],
                            style={
                                "display": "inline-block",
                                "width": "49%",
                                #'display': 'none',
                            },
                        ),
                        html.Div(
                            [dcc.Graph(id="metagene-plot")],
                            style={
                                "display": "inline-block",
                                "width": "49%",
                                #'display': 'none',
                            },
                        ),
                    ]
                )
            ]
        ),
    ]
)


@app.callback(
    Output("srametadata-table", "children"),
    [Input("srp", "value"), Input("assembly", "value")],
)
def display_datatable(srp, assembly):  # , gene_name):
    if isinstance(srp, dict):
        srp = srp["value"]
    df = get_srp_table(__DATASETS__, srp)
    return generate_table(df)  # .to_dict("rows")


@app.callback(Output("srp", "options"), [Input("assembly", "value")])
def update_project_dropdown(assembly):
    if isinstance(assembly, dict):
        assembly = assembly["value"]
    projects = get_projects(__DATASETS__, assembly)
    return projects


@app.callback(Output("srp", "value"), [Input("srp", "options")])
def update_srp_value(srp):
    if len(srp) != 0:
        return srp[0]
    else:
        return None


@app.callback(
    Output("read_length", "options"),
    [Input("srp", "value"), Input("assembly", "value")],
)
def generate_read_length_dropdown(srp, assembly):
    if isinstance(assembly, dict):
        assembly = assembly["value"]
    if isinstance(srp, dict):
        srp = srp["value"]
    return get_srp_read_lengths(__DATASETS__, srp)


@app.callback(
    Output("metagene-plot", "figure"),
    [
        Input("srp", "value"),
        Input("assembly", "value"),
        Input("normalize-chk", "values"),
        Input("read_length", "value"),
    ],
)
def display_metagene_plot(srp, assembly, normalize, length):  # gene_name, n_clicks):
    if isinstance(srp, dict):
        srp = srp["value"]
    metadata_tsv = get_project_summary_file(__DATASETS__, srp)
    metagene_dfs = project_summary_metagene_creator(metadata_tsv)
    return plot_metagene_coverage(metagene_dfs, length)


"""
@app.callback(
    Output("length-dist-plot", "figure"),
    [Input("srp", "value"), Input("assembly", "value")],
    # [State("gene-input", "value")],
    # [Input("gene-submit", "n_clicks")],
)
def display_read_length_dist_plot(srp, assembly):  # , state, n_clicks):
    if isinstance(srp, dict):
        srp = srp["value"]
    read_lengths = collect_read_lengths(srp, assembly)
    try:
        srp_df = get_srp_table_minimal(srp, assembly, "/data1")
        srp_df = (
            srp_df[["experiment_accession", "experiment_title"]]
            .drop_duplicates()
            .set_index("experiment_accession")
        )
        new_read_lengths = {
            key + "-" + str(srp_df.loc[key, "experiment_title"]): v
            for key, v in read_lengths.items()
        }
    except:
        new_read_lengths = read_lengths
    if not new_read_lengths:
        return
    return plot_read_length_distribution(
        new_read_lengths, shared_yaxes=False, plot_ridge=False
    )


@app.callback(
    Output("coherence-heatmap", "figure"),
    [Input("srp", "value"), Input("assembly", "value")],
    # [State("gene-input", "value")],
    # [Input("gene-submit", "n_clicks")],
)
def display_coherence_plot(srp, assembly):  # , state, n_clicks):
    if isinstance(srp, dict):
        srp = srp["value"]
    read_lengths = get_available_read_lengths(srp, assembly)
    df = get_coherence_scores(srp, assembly, read_lengths)
    return plot_coherence_heatmap(df)


"""
