import dash
import dash_core_components as dcc
import dash_html_components as html

# import dash_table_experiments as dt

from dash.dependencies import Output
from dash.dependencies import Input
from dash.dependencies import State

external_css = [
    "https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css",
    "https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
    "//fonts.googleapis.com/css?family=Raleway:400,300,600",
    "//fonts.googleapis.com/css?family=Dosis:Medium",
    "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
    "https://codepen.io/chriddyp/pen/bWLwgP.css",
    "https://codepen.io/chriddyp/pen/brPBPO.css",
    "https://unpkg.com/tachyons@4.10.0/css/tachyons.min.css",
]

app = dash.Dash(name="ribopod")  # title="ribopod : Visualizing ribo-seq datasets")
app.title = "ribopod : Visualizing ribo-seq datasets"

for css in external_css:
    app.css.append_css({"external_url": css})
# app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})
# Loading screen CSS
# app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/brPBPO.css"})
# Datatable CSS
# app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})

# Dash CSS
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})
# Loading screen CSS
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/brPBPO.css"})

external_js = ["https://code.jquery.com/jquery-3.2.1.min.js"]

for js in external_js:
    app.scripts.append_script({"external_url": js})

# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server
# app.scripts.config.serve_locally = True
app.config.suppress_callback_exceptions = True
