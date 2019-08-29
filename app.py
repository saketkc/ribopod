import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from app_multipage import app
from apps import landing_page, visualize_page
server = app.server

app.layout = html.Div(
    [dcc.Location(id="url", refresh=False), html.Div(id="page-content")]
)

external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]


@app.callback(Output("page-content", "children"), [Input("url", "pathname")])
def display_page(pathname):
    if pathname == "/":
        return landing_page.layout
    elif pathname == "/visualize":
        return visualize_page.layout
    else:
        return "404"


if __name__ == "__main__":
    app.run_server(host="68.181.32.41", debug=True, port=8050)
