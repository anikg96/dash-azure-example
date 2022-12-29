from dash import Dash, html, dcc, Input, Output, dash_table
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
import scanpy as sc

### Data Preprocessing ###
# Read the local.h5ad file from the 'Data' directory
adata = sc.read_h5ad('Data/local.h5ad')

# Create a dataframe from the UMAP parameters
umap_df = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"])

# Add the following columns to the UMAP dataframe
umap_df["cell_type"] = adata.obs["cell_type_original"].values
umap_df["cluster"] = adata.obs["final_cluster"].values
umap_df["disease"] = adata.obs["disease"].values

# The following two lines are very important.
# They are needed for changing the feature names from ENSG000___ to normal names
adata.var.reset_index(inplace=True)
adata.var.set_index('feature_name', inplace=True)

# Get the feature names
feature_names = adata.var.index.tolist()

# Create a data frame from the adata.X sparse matrix and label columns with feature names
# We set the adata.X type to bool to reduce its size
x_dataframe = pd.DataFrame.sparse.from_spmatrix(adata.X.astype(bool), columns=feature_names)

# Finally join the sparse matrix columns with the initial dataframe containing UMAP and cell_type information
final_data = pd.concat([umap_df, x_dataframe], axis=1)

# Create another dataframe from adata.X for the violin plots
# This time, we keep the original float32 datatype for the violin plots to show
violin_df = pd.DataFrame.sparse.from_spmatrix(adata.X, columns=feature_names)

### End of Data Preprocessing ###

app = Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    suppress_callback_exceptions=True
    # Might consider adding meta_tags later
)  # Initializes the app

### Creation of app components ###

## Common Layout Components ##
# Navbar
navbar = dbc.NavbarSimple(
    children=[
        dbc.NavItem(
            dbc.Button(
                "About",
                href="https://www.ukaachen.de/kliniken-institute/institut-fuer-experimentelle-innere-medizin-und-systembiologie/hayat-lab-in-translational-data-science/",
                color="primary",
                outline=True
            )
        )
    ],
    brand="Hayat Lab",
    brand_href="/",
    color="dark",
    dark=True
)

# Sidebar component and its style dict
SIDEBAR_STYLE = {
    "position": "fixed",
    "top": "7.8vh",
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#d7f5f5",
}
sidebar = html.Div(
    [
        html.H2("SciViewer", className="display-4"),
        html.Hr(),
        html.P(
            "Description for SciViewer App", className="lead"
        ),
        dbc.Nav(
            [
                dbc.NavLink("Page 1 - UMAP & Violin", href="/", active="exact"),
                dbc.NavLink("Page 2 - Table Select", href="/tableselect", active="exact"),
            ],
            vertical=True,
            pills=True,
        ),
        html.Footer("Copyright © Hayat Lab 2022")
    ],
    style=SIDEBAR_STYLE,
)

# Content component and its style dict
CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}
content = html.Div(id="page-content", style=CONTENT_STYLE)

# Combine them both in a div and call it sidebar and Content
# This step is done in the app.layout assignment itself

## End of Common Layout Components ##

## Page 1 ##
# Components: Main Graph, Reference Graph, Radio, Filter and Gene Dropdowns
# Main Graph
main_graph_figure = px.scatter(final_data, x="UMAP_1", y="UMAP_2", color="cell_type")
main_graph_component = dcc.Graph(figure=main_graph_figure, id="main-graph")

# Violin Plot
violin_plot = px.violin(violin_df, y="SNRPG")
violin_graph_component = dcc.Graph(figure=violin_plot, id="violin-plot")

# Radio group for three filtering options
radioItems = dcc.RadioItems([
    {'label': 'Cell Type', 'value': 'cell_type'},
    {'label': 'Disease', 'value': 'disease'},
    {'label': 'Cluster', 'value': 'cluster'}
],
    value='cell_type',
    id='filter-radio-page1'
)

# Dropdown of values for filtering based on selected radio button
filterDropdown = dcc.Dropdown(id='filter-dropdown-page1', clearable=True, searchable=True, multi=True)

# Multi-select dropdown for all genes
geneDropdown = dcc.Dropdown(
    options=[{'label': feature, 'value': feature} for feature in adata.var.index.unique()],
    clearable=True,
    searchable=True,
    id='gene-dropdown-page1',
    multi=True
)

# Final layout of Page 1 #
page1 = dbc.Container(
    [
        dbc.Row(
            dbc.Col(
                [
                    html.Label("Please use the following controls:"),
                    radioItems
                ],
                width=6
            ),
            justify="start"
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Label("Filter based on radio input:"),
                        filterDropdown
                    ],
                    width=6
                ),
                dbc.Col(
                    [
                        html.Label("Gene selection:"),
                        geneDropdown
                    ],
                    width=6
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col([
                    dcc.Loading(main_graph_component)
                ], width=6),
                dbc.Col([
                    dcc.Loading(violin_graph_component)
                ], width=6)
            ]
        )
    ],
    fluid=True
)
## End of Page 1 ##

## Page 2 ##
# UMAP
page2umap_figure = px.scatter(final_data[["UMAP_1", "UMAP_2", "cell_type"]], x="UMAP_1", y="UMAP_2", color="cell_type")
page2umap_component = dcc.Graph(figure=page2umap_figure, id="page2-umap")

# Violin Plot
page2_violin = px.violin(violin_df, y="SNRPG")
page2_violin_component = dcc.Graph(figure=page2_violin, id="page2-violin")

# Data Table
csv_table = pd.read_excel('aggregated_select_DE_genes.xlsx')
data_table = dash_table.DataTable(
    csv_table.to_dict('records'),
    [{"name": i.capitalize(), "id": i} for i in csv_table.columns],
    row_selectable="multi",
    fixed_columns={"headers": True},
    style_data={
        'whiteSpace': 'normal',
        'height': 'auto',
    },
    style_header={
        'backgroundColor': 'rgb(30, 30, 30)',
        'color': 'white'
    },
    style_table={'minWidth': '100%', 'overflowX': 'auto'},
    filter_action="native",
    sort_action="native",
    sort_mode="multi",
    cell_selectable=False,
    page_size=10,
    id="page2-datatable",
    #style_table={'overflowX': 'auto'}
)

# Final Page2 Layout
page2 = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col([
                    dcc.Loading(page2umap_component)
                ], width=6),
                dbc.Col([
                    dcc.Loading(page2_violin_component)
                ], width=6)
            ]
        ),
        dbc.Row(
            dbc.Col(data_table, width=12, style={"border-style": "solid"}),
        ),
        dbc.Alert(id='absent-genes')
    ],
    fluid=True
)
## End of Page 2 ##
### Layout for the app ###
app.layout = html.Div(
    [
        dcc.Location(id="url"),
        navbar,
        sidebar,
        dcc.Loading(content, id="content-loading", type="circle")
    ]
)
# Needs serious work
# app.layout = html.Div([
#     navbar,
#     radioItems,
#     filterDropdown,
#     dcc.Loading(main_graph_component, id="main-graph-loading", type="default"),
#     geneDropdown
# ])
### End of Layout ###

### Callbacks begin here ###

# Call back for navigating between two pages
@app.callback(
    Output("page-content", "children"),
    Input("url", "pathname")
)
def render_page_content(pathname):
    if pathname == "/":
        return page1
    elif pathname == "/tableselect":
        return page2
    # If the user tries to reach a different page, return a 404 message
    return html.Div(
        [
            html.H1("404: Not found", className="text-danger"),
            html.Hr(),
            html.P(f"The pathname {pathname} was not recognised..."),
        ],
        className="p-3 bg-light rounded-3",
    )

# Callback for populating filter dropdown based on selected radio button
@app.callback(
    Output('filter-dropdown-page1', 'options'),
    Output('filter-dropdown-page1', 'value'),
    Input('filter-radio-page1', 'value')
)
def populate_dropdown_page1(radiovalue):
    dropdown_list = [{'label': value, 'value': value} for value in umap_df[radiovalue].values.unique()]
    return dropdown_list, None # Setting the dropdown's value to 'None' fixes a KeyError Bug

# Callback for filtering the graph based on (cell-type/disease/final cluster) filter
@app.callback(
    Output('main-graph', 'figure'),
    Input('filter-radio-page1', 'value'),
    Input('filter-dropdown-page1', 'value'),
    Input('gene-dropdown-page1', 'value'),
)
def update_page1_main_graph_from_controls(radiovalue, filtervalue, genevalues):
    if not filtervalue:
        if not genevalues:
            filtered_data1 = final_data[["UMAP_1", "UMAP_2", radiovalue]]
            figure1 = px.scatter(data_frame=filtered_data1, x="UMAP_1", y="UMAP_2", color=radiovalue)
            return figure1
        elif genevalues:
            filtered_data2 = final_data[genevalues + ["UMAP_1", "UMAP_2", radiovalue]]
            for gene in genevalues:
                mask = filtered_data2[gene] == 1
                filtered_data2 = filtered_data2[mask]
            figure2 = px.scatter(data_frame=filtered_data2, x="UMAP_1", y="UMAP_2", color=filtered_data2[radiovalue].values.astype(str))
            return figure2
    elif not genevalues:
        filtered_data3 = final_data[["UMAP_1", "UMAP_2", radiovalue]]
        filtered_data3 = filtered_data3[filtered_data3[radiovalue].isin(filtervalue)]
        figure3 = px.scatter(data_frame=filtered_data3, x="UMAP_1", y="UMAP_2", color=filtered_data3[radiovalue].values.astype(str))
        return figure3
    else:
        filtered_data4 = final_data[genevalues + ["UMAP_1", "UMAP_2", radiovalue]]
        filtered_data4 = filtered_data4[filtered_data4[radiovalue].isin(filtervalue)]
        for gene in genevalues:
            mask = filtered_data4[gene]==1
            filtered_data4 = filtered_data4[mask]
        figure4 = px.scatter(data_frame=filtered_data4, x="UMAP_1", y="UMAP_2", color=filtered_data4[radiovalue].values.astype(str))
        return figure4

# Update Violin Plot
@app.callback(
    Output("violin-plot", "figure"),
    Input("gene-dropdown-page1", "value")
)
def update_violin_plot(genevalues):
    if not genevalues:
        return px.violin(violin_df[["CSTF2T", "ALDH1A2", "SNRPG"]], y=["CSTF2T", "ALDH1A2", "SNRPG"])
    else:
        return px.violin(violin_df[genevalues], y=genevalues)

# Page 2 Callback - Update both plots when a row is selected
@app.callback(
    Output("page2-umap", "figure"),
    Output("page2-violin", "figure"),
    Output("absent-genes", "children"),
    Input("page2-datatable", "selected_rows")
)
def update_page2_plots(selected_rows):
    if selected_rows:
        selected_genes = csv_table.iloc[selected_rows]["gene"]
        valid_genes = [gene for gene in selected_genes if gene in feature_names]
        valid_genes = set(valid_genes)
        absent_genes = list(set(selected_genes) - valid_genes)
        valid_genes = list(valid_genes)
        filtered_data = final_data[valid_genes + ["UMAP_1", "UMAP_2", "cell_type"]]
        for gene in valid_genes:
            mask = filtered_data[gene] == 1
            filtered_data = filtered_data[mask]
        umap_updated = px.scatter(data_frame=filtered_data, x="UMAP_1", y="UMAP_2", color=filtered_data.values.astype(str))
        violin_updated = px.violin(violin_df[valid_genes], y=valid_genes)
        return umap_updated, violin_updated, "Following genes missing: " + ','.join(absent_genes)
    else:
        umap_original = px.scatter(final_data[["UMAP_1", "UMAP_2", "cell_type"]], x="UMAP_1", y="UMAP_2", color="cell_type")
        violin_original = px.violin(violin_df, y="SNRPG")
        return umap_original, violin_original, "No genes selected"

### Callbacks end here ###

# Run the app on port 8050. Remove debug=True for prod
if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
