import json
from app_multipage import app

import dash_core_components as dcc
import dash_html_components as html

from dash.dependencies import Input, Output
import flask
from flask import send_file
from plotly import tools
from plotly.graph_objs import Bar

from init import __DATASETS__
from init import __SPECIES__
from init import build_to_species

import pandas as pd
from plotly.graph_objs import Bar, Figure, Layout

datasets_counts = __DATASETS__.index.value_counts().to_dict()
print(datasets_counts)
datasets_counts = {
    build_to_species(key): value for key, value in datasets_counts.items()
}
datasets_counts = pd.Series(datasets_counts).sort_values()


def display_bar_plot():  # , state, n_clicks):
    trace = Bar(
        x=datasets_counts.index.tolist(),
        y=datasets_counts.values.tolist(),
        name="Number of Datasets",
    )
    # layout = Layout(clickmode='event+select'),
    fig = Figure(data=[trace])  #
    fig["layout"].update(clickmode="event+select")
    # print(fig)
    return fig  # layout=layout)


figure = display_bar_plot()
layout = html.Div(
    [
        html.Div(
            [
                html.Div(
                    [
                        html.H2(
                            "ribopod: A database of actively-translating ORFs in Ribo-seq data"
                        )
                    ],
                    style={
                        "margin-right": "auto",
                        "margin-left": "auto",
                        "width": "100%",
                    },
                ),
                html.Hr(),
                html.Div(
                    [
                        html.Div(
                            [html.A("      Available Datasets", href="/")],
                            style={"display": "inline-block", "width": "30%"},
                        ),
                        html.Div(
                            [html.A("      Visualize", href="/visualize")],
                            style={"display": "inline-block", "width": "30%"},
                        ),
                    ]
                ),
                html.Hr(),
            ]
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                dcc.Loading(
                                    id="loading-dataset-distribution",
                                    children=[
                                        html.Div(
                                            [
                                                dcc.Graph(
                                                    id="dataset-distribution",
                                                    figure=figure,
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
                ),
                html.Div(id="species-data"),
            ]
        ),
    ]
)


@app.callback(
    Output("species-data", "children"), [Input("dataset-distribution", "selectedData")]
)
def display_click_data(click_data):
    data = json.dumps(click_data, indent=2)
    print(data)
    return data
