import dash_core_components as dcc
import dash_html_components as html

from dash.dependencies import Input, Output
import flask
from flask import send_file

from init import __DATASETS__
from init import __SPECIES__

from app_multipage import app
from dash_helper import generate_table, path_leaf
from fragment_length_helper import (
    parse_ribotricer_bam_summary,
    project_summary_read_length_creator,
    plot_read_length_distribution,
)
from metagene_helper import (
    project_summary_metagene_creator,
    plot_metagene_coverage,
    metagene_profile_to_phase_score_matrix,
    plot_phase_score_heatmap,
)

from project_helper import (
    get_projects,
    get_srp_table,
    get_srp_read_lengths,
    get_project_summary_file,
    get_summarized_phase_scores,
    get_summarized_orf_counts,
)

from orf_helper import plot_orf_counts_stacked_bar, plot_phase_scores_violin

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
        html.Div(id="datatable-interactivity-container"),
        html.Div(
            [
                html.Div(
                    [
                        html.Label("Species: "),
                        dcc.Dropdown(
                            options=__SPECIES__,
                            value="SC5314",
                            id="assembly",
                            searchable=True,
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
                        html.Label("SRP: "),
                        dcc.Dropdown(options=[], value=[], id="srp", searchable=True),
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
                                html.Label("Fragment length for metagene plot: "),
                                dcc.Dropdown(
                                    options=[],
                                    value=28,
                                    searchable=True,
                                    id="read_length",
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
                                            "value": "True",
                                        }
                                    ],
                                    id="normalize-chk",
                                    value=[],
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
                        dcc.Tabs(
                            id="tabs-div",
                            value="tab-phase-score",
                            children=[
                                dcc.Tab(
                                    label="Phase Score",
                                    value="tab-phase-score",
                                    children=[
                                        html.Div(
                                            [
                                                html.Div(
                                                    [
                                                        dcc.Loading(
                                                            id="loading-coherence-heatmap",
                                                            children=[
                                                                html.Div(
                                                                    [
                                                                        dcc.Graph(
                                                                            id="coherence-heatmap",
                                                                            style={
                                                                                "margin-right": "auto",
                                                                                "margin-left": "auto",
                                                                                "width": "100%",
                                                                            },
                                                                        )
                                                                    ]
                                                                )
                                                            ],
                                                            type="circle",
                                                        )
                                                    ],
                                                    style={"width": "100%"},
                                                )
                                            ]
                                        )
                                    ],
                                ),
                                dcc.Tab(
                                    label="Distribution of ORFs",
                                    value="tab-orf-dist",
                                    children=[
                                        html.Div(
                                            [
                                                html.Div(
                                                    [
                                                        dcc.Loading(
                                                            id="loading-orf-count-dist-plot",
                                                            children=[
                                                                html.Div(
                                                                    [
                                                                        dcc.Graph(
                                                                            id="orf-count-dist-plot",
                                                                            style={
                                                                                "margin-right": "auto",
                                                                                "margin-left": "auto",
                                                                                "width": "100%",
                                                                            },
                                                                        )
                                                                    ]
                                                                )
                                                            ],
                                                            type="circle",
                                                        )
                                                    ],
                                                    style={"width": "100%"},
                                                )
                                            ]
                                        )
                                    ],
                                ),
                                dcc.Tab(
                                    label="Distribution of Phase Scores",
                                    value="tab-phase-score-dist",
                                    children=[
                                        html.Div(
                                            [
                                                html.Div(
                                                    [
                                                        dcc.Loading(
                                                            id="loading-phase-scire-dist-plot",
                                                            children=[
                                                                html.Div(
                                                                    [
                                                                        dcc.Graph(
                                                                            id="phase-score-dist-plot",
                                                                            style={
                                                                                "margin-right": "auto",
                                                                                "margin-left": "auto",
                                                                                "width": "100%",
                                                                            },
                                                                        )
                                                                    ]
                                                                )
                                                            ],
                                                            type="circle",
                                                        )
                                                    ],
                                                    style={"width": "100%"},
                                                )
                                            ]
                                        )
                                    ],
                                ),
                                dcc.Tab(
                                    label="Metadata",
                                    value="tab-sra-metadata",
                                    children=[
                                        html.Div(
                                            html.Div(
                                                id="srametadata-table",
                                                style={
                                                    # "overflowX": "scroll",
                                                    # "overflowY": "scroll",
                                                    "margin-right": "auto",
                                                    "margin-left": "auto",
                                                    "width": "100%",
                                                    "text-align": "center",
                                                },
                                            ),
                                            style={
                                                "width": "100%",
                                                "height": "30%",
                                                "display": "inline-block",
                                            },
                                        )
                                    ],
                                ),
                            ],
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
                            [
                                dcc.Loading(
                                    id="loading-length-dist-plot",
                                    children=[
                                        html.Div([dcc.Graph(id="length-dist-plot")])
                                    ],
                                    type="graph",
                                )
                            ],
                            style={
                                "display": "inline-block",
                                "width": "49%",
                                #'display': 'none',
                            },
                        ),
                        html.Div(
                            [
                                dcc.Loading(
                                    id="loading-metagene-plot",
                                    children=[
                                        html.Div([dcc.Graph(id="metagene-plot")])
                                    ],
                                    type="graph",
                                )
                            ],
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
        return srp[0]["value"]
    else:
        return []


@app.callback(
    Output("read_length", "options"),
    [Input("srp", "value"), Input("assembly", "value")],
)
def generate_read_length_dropdown(srp, assembly):
    if isinstance(assembly, dict):
        assembly = assembly["value"]
    if isinstance(srp, dict):
        srp = srp["value"]
    options = get_srp_read_lengths(__DATASETS__, srp)
    return options


@app.callback(
    Output("metagene-plot", "figure"),
    [
        Input("srp", "value"),
        Input("assembly", "value"),
        Input("normalize-chk", "value"),
        Input("read_length", "value"),
    ],
)
def display_metagene_plot(srp, assembly, normalize, length):  # gene_name, n_clicks):
    if isinstance(srp, dict):
        srp = srp["value"]
    project_summary_file = get_project_summary_file(__DATASETS__, srp, assembly)
    metagene_dfs = project_summary_metagene_creator(project_summary_file)
    return plot_metagene_coverage(metagene_dfs, length, normalize_per_codon=normalize)


@app.callback(
    Output("length-dist-plot", "figure"),
    [Input("srp", "value"), Input("assembly", "value")],
)
def display_read_length_dist_plot(srp, assembly):
    if isinstance(srp, dict):
        srp = srp["value"]
    if isinstance(assembly, dict):
        assembly = assembly["value"]
    project_summary_file = get_project_summary_file(__DATASETS__, srp, assembly)
    read_length_dist = project_summary_read_length_creator(project_summary_file)
    return plot_read_length_distribution(
        read_length_dist, shared_yaxes=False, plot_ridge=False
    )


@app.server.route("/download")
def serve_static():
    path = flask.request.args.get("value")
    filename = path_leaf(path)
    return send_file(
        path, mimetype="text/csv", attachment_filename=filename, as_attachment=True
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
    if isinstance(assembly, dict):
        assembly = assembly["value"]
    project_summary_file = get_project_summary_file(__DATASETS__, srp, assembly)
    metagene_dfs = project_summary_metagene_creator(project_summary_file)
    phase_score_df = metagene_profile_to_phase_score_matrix(metagene_dfs)
    return plot_phase_score_heatmap(phase_score_df)


@app.callback(
    Output("orf-count-dist-plot", "figure"),
    [Input("srp", "value"), Input("assembly", "value")],
)
def display_orf_count_dist_plot(srp, assembly):  # , state, n_clicks):
    if isinstance(srp, dict):
        srp = srp["value"]
    if isinstance(assembly, dict):
        assembly = assembly["value"]
    # get_summarized_phase_scores,
    orf_df = get_summarized_orf_counts(__DATASETS__, srp, assembly)
    if orf_df:
        return plot_orf_counts_stacked_bar(orf_df)


@app.callback(
    Output("phase-score-dist-plot", "figure"),
    [Input("srp", "value"), Input("assembly", "value")],
)
def display_phase_score_dist_plot(srp, assembly):  # , state, n_clicks):
    if isinstance(srp, dict):
        srp = srp["value"]
    if isinstance(assembly, dict):
        assembly = assembly["value"]
    # get_summarized_phase_scores,
    phase_scores_df = get_summarized_phase_scores(__DATASETS__, srp, assembly)
    return plot_phase_scores_violin(phase_scores_df)
