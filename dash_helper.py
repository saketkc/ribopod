import ntpath
import dash_html_components as html


def path_leaf(path):
    """Get path's tail from a filepath"""
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def generate_table(dataframe):
    """Generate dash table for pandas dataframe.

    Parameters
    ----------
    dataframe: pandas.DataFrame

    Returns
    -------
    table: html.Table
    """
    rows = []
    cols_to_retain = [
        "experiment_accession",
        "experiment_title",
        "run_accession",
        "ribotricer_orf",
        "ribotricer_metagene_5p",
        "ribotricer_metagene_3p",
        "ribotricer_metagene_plot",
        "ribotricer_protocol",
        "ribotricer_bam_summary",
    ]
    dataframe = dataframe.loc[:, cols_to_retain]
    for i in range(len(dataframe)):
        row = []
        for col in dataframe.columns:
            value = dataframe.iloc[i][col]
            if (
                col
                in [
                    "ribotricer_orfs",
                    "ribotricer_metagene_5p",
                    "ribotricer_metagene_3p",
                    "ribotricer_metagene_plot",
                    "ribotricer_protocol",
                    "ribotricer_bam_summary",
                ]
                and str(value) != "nan"
            ):
                try:
                    cell = html.Td(
                        html.A(
                            href="/download?value=" + value,
                            children=[html.I(className="fa fa-download")],
                        ),
                        style={"text-align": "center"},
                    )
                except:
                    print(value)
            else:
                cell = html.Td(children=value)
            row.append(cell)
        rows.append(html.Tr(row))
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])] + rows,
        style={
            "width": "100%",
            "margin-right": "auto",
            "margin-left": "25%",
            "align": "center",
            "display": "inline-block",
        },
    )
