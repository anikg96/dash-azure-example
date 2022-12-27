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

### Layout for the app ###
# Needs serious work
app.layout = html.Div([
    radioItems,
    filterDropdown,
    dcc.Loading(main_graph_component, id="main-graph-loading", type="default"),
    geneDropdown,
    #html.Div(id="sample-div", style={'border-style': 'solid'}),
    #html.Div(id="sample-div2", style={'border-style': 'solid'})
])
### End of Layout ###

### Callbacks begin here ###

# Callback for populating filter dropdown based on selected radio button
@app.callback(
    Output('filter-dropdown-page1', 'options'),
    #Output('filter-dropdown-page1', 'value'),
    Input('filter-radio-page1', 'value')
)
def populate_dropdown_page1(radiovalue):
    dropdown_list = [{'label': value, 'value': value} for value in umap_df[radiovalue].values.unique()]
    #dropdown_list.insert(0, {'label': 'Please select a value', 'value': 'Nothing'})
    return dropdown_list#, 'Nothing'

# Callback for filtering the graph based on (cell-type/disease/final cluster) filter
@app.callback(
    Output('main-graph', 'figure'),
    #Output('sample-div', 'children'),
    #Output('sample-div2', 'children'),
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
        figure3 = px.scatter(data_frame=filtered_data3, x="UMAP_1", y="UMAP_2", color=filtered_data3["cell_type"].values.astype(str))
        return figure3
    else:
        filtered_data4 = final_data[genevalues + ["UMAP_1", "UMAP_2", radiovalue]]
        filtered_data4 = filtered_data4[filtered_data4[radiovalue].isin(filtervalue)]
        for gene in genevalues:
            mask = filtered_data4[gene]==1
            filtered_data4 = filtered_data4[mask]
        figure4 = px.scatter(data_frame=filtered_data4, x="UMAP_1", y="UMAP_2", color=filtered_data4["cell_type"].values.astype(str))
        # Remove these maybe?
        return figure4

### Callbacks end here ###

if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
