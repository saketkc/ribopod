import dash_html_components as html


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
    for i in range(len(dataframe)):
        row = []
        for col in dataframe.columns:
            value = dataframe.iloc[i][col]
            # update this depending on which
            # columns you want to show links for
            # and what you want those links to be
            if col == "ribotricer_orfs" and value:
                cell = html.Td(html.A(href=value, children="Download"))
            else:
                cell = html.Td(children=value)
            row.append(cell)
        rows.append(html.Tr(row))
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])]
        + rows
    )
