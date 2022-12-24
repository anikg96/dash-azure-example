from dash import Dash, html, dcc, Input, Output
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

app = Dash(__name__)  # Initializes the app

### Creation of app components ###

# Main Graph
main_graph_figure = px.scatter(final_data, x="UMAP_1", y="UMAP_2", color="cell_type")
main_graph_component = dcc.Graph(figure=main_graph_figure, id="main-graph")

# Radio group for three filtering options
radioItems = dcc.RadioItems([
    {'label': 'Cell Type', 'value': 'cell_type'},
    {'label': 'Disease', 'value': 'disease'},
    {'label': 'Cluster', 'value': 'cluster'}
],
    value='cell_type',
    id='filter_radio'
)

# Dropdown of values for filtering based on selected radio button
filterDropdown = dcc.Dropdown(id='filter-dropdown-page1', clearable=False, searchable=True)

# Multi-select dropdown for all genes
geneDropdown = dcc.Dropdown(
    options=[{'label': feature, 'value': feature} for feature in adata.var.index.unique()],
    clearable=True,
    searchable=True,
    id='gene-dropdown-page1',
    multi=True
)

### Layout for the app ###
# Needs serious work
app.layout = html.Div([
    radioItems,
    filterDropdown,
    main_graph_component,
    geneDropdown
])
### End of Layout ###

### Callbacks begin here ###

# Callback for populating filter dropdown based on selected radio button
@app.callback(
    Output('filter-dropdown-page1', 'options'),
    Output('filter-dropdown-page1', 'value'),
    Input('filter_radio', 'value')
)
def populate_dropdown_page1(radiovalue):
    dropdown_list = [{'label': value, 'value': value} for value in umap_df[radiovalue].values.unique()]
    dropdown_list.insert(0, {'label': 'Please select a value', 'value': 'Nothing'})
    return dropdown_list, 'Nothing'

# Callback for filtering the graph based on (cell-type/disease/final cluster) filter
# @app.callback(
#     Output('main-graph', 'figure'),
#
# ) Work on this later
### Callbacks end here ###

if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
