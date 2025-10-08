# app.py
import dash
from dash import dcc, html, dash_table, Dash
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import sqlite3
import numpy as np
import os
from utils import *

CACHE_DIR = 'cache'
DATA_DIR = 'data/'
os.makedirs(CACHE_DIR, exist_ok=True)

# Read in data for targeted screen
cell_lines = pd.DataFrame({
    "individual": ["eipl", "eipl", "iudw", "iudw", "jejf", "jejf", "kolf", "kolf", "paab", "paab"],
    "Cell_Line": ["eipl_1", "eipl_3", "iudw_1", "iudw_4", "jejf_2", "jejf_3", "kolf_2", "kolf_3", "paab_3", "paab_4"],
})
targeted_screen = query_database("SELECT * FROM targeted_screen", CACHE_DIR, DATA_DIR)
targeted_screen = pd.merge(targeted_screen, cell_lines, on="Cell_Line", how="left")
targeted_screen['linetype'] = np.where(targeted_screen['pval_adj'] < 0.1, "solid", "dash")
heritability = pd.read_csv(DATA_DIR + "ST3-3_Heritable_Trans_Effects.tsv", sep="\t")

# Using a static color scale 
colors = px.colors.diverging.RdBu[::-1]
n_colors = 101
color_scale = [colors[int(i * (len(colors) - 1) / (n_colors - 1))] for i in range(n_colors)]

# Get unique targets and genes per data frame
target_names_lfc = query_database("SELECT DISTINCT Target FROM lfc", CACHE_DIR, DATA_DIR)['Target'].tolist()
#target_names_lfc = pd.read_sql_query("SELECT DISTINCT Target FROM lfc", con)['Target'].tolist()
gene_names_lfc = query_database("SELECT DISTINCT Expressed_Gene_Symbol FROM lfc", CACHE_DIR, DATA_DIR)['Expressed_Gene_Symbol'].tolist()
target_names_cor = query_database("SELECT DISTINCT Target_1 FROM target_cor", CACHE_DIR, DATA_DIR)['Target_1'].tolist()
gene_names_cor = query_database("SELECT DISTINCT Expressed_Gene_1 FROM gene_cor", CACHE_DIR, DATA_DIR)['Expressed_Gene_1'].tolist()
target_names_tar = query_database("SELECT DISTINCT Target FROM targeted_screen", CACHE_DIR, DATA_DIR)['Target'].tolist()
gene_names_tar = query_database("SELECT DISTINCT Expressed_Gene_Symbol FROM targeted_screen", 
                                CACHE_DIR, DATA_DIR)['Expressed_Gene_Symbol'].tolist()
target_names_heritability = heritability['Target'].unique().tolist()
gene_names_heritability = heritability['Downstream_Gene'].unique().tolist()

# Define UI for app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.YETI])

app.layout = dbc.Container(
    [
        dcc.Tabs(
            id="tabs",
            children=[
                # Tab 1: LFC Heatmap
                dcc.Tab(
                    label="Log fold change",
                    children=[
                        dbc.Card(
                            [
                                dbc.CardHeader("Log fold change heatmap"),
                                dbc.CardBody(
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                [
                                                    html.P("Select targets:"),
                                                    dcc.Dropdown(
                                                        id="selected_targets",
                                                        options=[{'label': i, 'value': i} for i in target_names_lfc],
                                                        value=["ARNT", "HIF1A", "VHL"],
                                                        multi=True,
                                                        placeholder='Type to search target genes...',
                                                    ),
                                                    html.Br(),
                                                    html.P("Select expressed genes:"),
                                                    dcc.Dropdown(
                                                        id="selected_genes",
                                                        options=[{'label': i, 'value': i} for i in gene_names_lfc],
                                                        value=["ALDOA", "BNIP3", "ENO1", "GAPDH", "LDHA", "PGAM1", "PGK1", "PKM"],
                                                        multi=True,
                                                        placeholder='Type to search expressed genes...',
                                                    ),
                                                    html.Br(),
                                                ],
                                                md=4,
                                            ),
                                            dbc.Col(
                                                [
                                                    html.Div(id="lfc_message"),
                                                    dcc.Graph(id="lfc_heatmap"),
                                                ],
                                                md=8,
                                            ),
                                        ]
                                    )
                                ),
                            ]
                        ),
                        dbc.Card(
                            [
                                dbc.CardHeader("Significant regulators for a selected expressed gene:"),
                                dbc.CardBody(
                                    [
                                        dbc.Row(
                                            [
                                                dbc.Col(
                                                    dcc.Dropdown(
                                                        id="selected_lfc_gene",
                                                        options=[{'label': i, 'value': i} for i in gene_names_lfc],
                                                        value="PGK1",
                                                        multi=False,
                                                        placeholder='Type to search expressed gene...',
                                                    ),
                                                    md=6,
                                                ),
                                                dbc.Col(
                                                    [
                                                        dcc.Slider(
                                                            id="num_regulators",
                                                            min=5,
                                                            max=50,
                                                            step=1,
                                                            value=10,
                                                            marks={i: str(i) for i in range(5, 51, 5)},
                                                        ),
                                                        html.Div(id="num_regulators_output", style={'fontSize': 12})
                                                    ],
                                                    md=6,
                                                ),
                                            ]
                                        ),
                                        html.Br(),
                                        dash_table.DataTable(id="table_sig_regulators")
                                    ],
                                ),
                            ],
                        ),
                        dbc.Card(
                            [
                                dbc.CardHeader("Significant genes for a selected perturbation:"),
                                dbc.CardBody(
                                    [
                                        dbc.Row(
                                            [
                                                dbc.Col(
                                                    dcc.Dropdown(
                                                        id="selected_lfc_target",
                                                        options=[{'label': i, 'value': i} for i in target_names_lfc],
                                                        value="HIF1A",
                                                        multi=False,
                                                        placeholder='Type to search target...',
                                                    ),
                                                    md=6,
                                                ),
                                                dbc.Col(
                                                    [
                                                        dcc.Slider(
                                                            id="num_genes",
                                                            min=5,
                                                            max=50,
                                                            step=1,
                                                            value=10,
                                                            marks={i: str(i) for i in range(5, 51, 5)},
                                                        ),
                                                        html.Div(id="num_genes_output", style={'fontSize': 12})
                                                    ],
                                                    md=6,
                                                ),
                                            ]
                                        ),
                                        html.Br(),
                                        dash_table.DataTable(id="table_sig_genes")
                                    ],
                                ),
                            ],
                        )
                    ]
                ),
                # Tab 2: Target-target correlation
                dcc.Tab(
                    label="Target-target correlation",
                    children=[
                        dbc.Card(
                            [
                                dbc.CardHeader("Target-target correlation heatmap"),
                                dbc.CardBody(
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                [
                                                    html.P("Select targets:"),
                                                    dcc.Dropdown(
                                                        id="selected_cor_targets",
                                                        options=[{'label': i, 'value': i} for i in target_names_cor],
                                                        value=["CTR9", "PAF1", "MED7", "MED9", "MED28", "SUPT20H", "CNOT1", "MAX", "ING3"],
                                                        multi=True,
                                                        placeholder='Type to search targets...',
                                                    )
                                                ],
                                                md=4,
                                            ),
                                            dbc.Col(
                                                [
                                                    html.Div(id="target_cor_message"),
                                                    dcc.Graph(id="target_cor_heatmap"),
                                                ],
                                                md=4,
                                            ),
                                        ]
                                    )
                                ),
                            ],
                        ),
                        dbc.Card(
                            [
                                dbc.CardHeader("Scatterplot for a selected target pair"),
                                dbc.CardBody(
                                    [
                                        dbc.Row(
                                            [
                                                dbc.Col(
                                                    [
                                                        html.P("Select target 1:"),
                                                        dcc.Dropdown(
                                                            id="selected_scatter_target1",
                                                            options=[{'label': i, 'value': i} for i in target_names_lfc],
                                                            value="CTR9",
                                                            multi=False,
                                                            placeholder='Type to search targets...',
                                                        ),
                                                    ],
                                                    md=3,
                                                ),
                                                dbc.Col(
                                                    [
                                                        html.P("Select target 2:"),
                                                        dcc.Dropdown(
                                                            id="selected_scatter_target2",
                                                            options=[{'label': i, 'value': i} for i in target_names_lfc],
                                                            value="PAF1",
                                                            multi=False,
                                                            placeholder='Type to search targets...',
                                                        ),
                                                    ],
                                                    md=3,
                                                ),
                                            ]
                                        ),
                                        html.Br(),
                                        dbc.Col(
                                            [
                                                dcc.Graph(id="target_scatter")
                                            ], 
                                            md = 6
                                        )
                                    ],
                                ),
                            ],
                        )
                    ]
                ),
                # Tab 3: Gene-gene correlation
                dcc.Tab(
                    label="Gene-gene correlation",
                    children=[
                        dbc.Card(
                            [
                                dbc.CardBody(
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                [
                                                    html.P("Select expressed genes:"),
                                                    dcc.Dropdown(
                                                        id="selected_cor_genes",
                                                        options=[{'label': i, 'value': i} for i in gene_names_cor],
                                                        value=["TERF1", "TDGF1", "CD24", "DNMT3B", "L1TD1"],
                                                        multi=True,
                                                        placeholder='Type to search targets...',
                                                    )
                                                ],
                                                md=4,
                                            ),
                                            dbc.Col(
                                                [
                                                    html.Div(id="gene_cor_message"),
                                                    dcc.Graph(id="gene_cor_heatmap"),
                                                ],
                                                md=4,
                                            ),
                                        ]
                                    )
                                ),
                            ],
                        ),
                        dbc.Card(
                            [
                                dbc.CardHeader("Scatterplot for a selected expressed gene pair"),
                                dbc.CardBody(
                                    [
                                        dbc.Row(
                                            [
                                                dbc.Col(
                                                    [
                                                        html.P("Select gene 1:"),
                                                        dcc.Dropdown(
                                                            id="selected_scatter_gene1",
                                                            options=[{'label': i, 'value': i} for i in gene_names_lfc],
                                                            value="TERF1",
                                                            multi=False,
                                                            placeholder='Type to search genes...',
                                                        )
                                                    ],
                                                    md=3,
                                                ),
                                                dbc.Col(
                                                    [
                                                        html.P("Select gene 2:"),
                                                        dcc.Dropdown(
                                                            id="selected_scatter_gene2",
                                                            options=[{'label': i, 'value': i} for i in gene_names_lfc],
                                                            value="CD24",
                                                            multi=False,
                                                            placeholder='Type to search genes...',
                                                        ),
                                                    ],
                                                    md=3,
                                                ),
                                            ]
                                        ),
                                        html.Br(),
                                        dbc.Col(
                                            [
                                                dcc.Graph(id="gene_scatter")
                                            ],
                                            md = 6
                                        )
                                    ],
                                ),
                            ],
                        )
                    ]
                ),
                # Tab 4: Targeted screen
                dcc.Tab(
                    label="Heritable effects",
                    children=[
                        dbc.Card(
                            [
                                dbc.CardHeader("Influences of genetic background on perturbation response"),
                                dbc.CardBody(
                                    [
                                        html.P("Visualization of cell-line specific effects for pairs with heritable effects (see tables below)"),
                                        dbc.Row(
                                            [
                                                dbc.Col(
                                                    [
                                                        html.P("Select target:"),
                                                        dcc.Dropdown(
                                                            id="selected_tar_target",
                                                            options=[{'label': i, 'value': i} for i in target_names_tar],
                                                            value="CREB3L2",
                                                            multi=False,
                                                            placeholder='Type to search expressed gene...',
                                                        ),
                                                    ],
                                                    md=3,
                                                ),
                                                dbc.Col(
                                                    [
                                                        html.P("Select expressed gene:"),
                                                        dcc.Dropdown(
                                                            id="selected_tar_gene",
                                                            options=[{'label': i, 'value': i} for i in gene_names_tar],
                                                            value="MICU2",
                                                            multi=False,
                                                            placeholder='Type to search expressed gene...',
                                                        ),
                                                    ],
                                                    md=3,
                                                ),
                                            ]
                                        ),
                                        html.Br(),
                                        dbc.Col(
                                            [
                                                dcc.Graph(id="cell_line_effects"),
                                            ],
                                            md = 6
                                        )
                                    ],
                                ),
                            ],
                        ),
                        dbc.Card(
                            [
                                dbc.CardHeader("Heritable trans effects for a selected expressed gene:"),
                                dbc.CardBody(
                                    [
                                        dbc.Row(
                                            dbc.Col(
                                                dcc.Dropdown(
                                                    id="selected_heritable_gene",
                                                    options=[{'label': i, 'value': i} for i in gene_names_heritability],
                                                    value="MICU2",
                                                    multi=False,
                                                    placeholder='Type to search expressed gene...',
                                                ),
                                                md=3,
                                            )
                                        ),
                                        html.Br(),
                                        dash_table.DataTable(id="table_heritable_gene")
                                    ],
                                ),
                            ],
                        ),
                        dbc.Card(
                            [
                                dbc.CardHeader("Heritable trans effects for a selected target (perturbation):"),
                                dbc.CardBody(
                                    [
                                        dbc.Row(
                                            dbc.Col(
                                                dcc.Dropdown(
                                                    id="selected_heritable_target",
                                                    options=[{'label': i, 'value': i} for i in target_names_heritability],
                                                    value="CREB3L2",
                                                    multi=False,
                                                    placeholder='Type to search target...',
                                                ),
                                                md=3,
                                            )
                                        ),
                                        html.Br(),
                                        dash_table.DataTable(id="table_heritable_target")
                                    ],
                                ),
                            ],
                        ),
                    ]
                ),
            ]
        )
    ],
    fluid=True,
)

# Define server logic
@app.callback(
    Output("lfc_heatmap", "figure"),
    Output("lfc_message", "children"),
    Input("selected_targets", "value"),
    Input("selected_genes", "value"),
)
def update_lfc_heatmap(selected_targets, selected_genes):
    if not selected_targets or not selected_genes:
        return go.Figure(), "Please select at least one target and one gene."
    if len(selected_targets) > 25 or len(selected_genes) > 25:
        return go.Figure(), "Only up to 25 targets and up to 25 expressed genes can be visualized at once. Please select 1-25 targets and 1-25 expressed genes for the heatmap."
    # Get LFC for selection
    filtered_lfc = get_filtered_lfc(selected_targets, selected_genes, CACHE_DIR, DATA_DIR)
    if filtered_lfc.empty:
        return go.Figure(), "No data available for the selected targets and genes."
    # Plot the heatmap for the selected data
    heatmap = create_lfc_heatmap(filtered_lfc, color_scale)
    return heatmap, ""

@app.callback(
    Output("table_sig_regulators", "data"),
    Input("selected_lfc_gene", "value"),
    Input("num_regulators", "value"),
)
def update_sig_regulators_table(selected_lfc_gene, num_regulators):
    query = f"""
        SELECT Target, lfc, pval_adj
        FROM lfc
        WHERE Expressed_Gene_Symbol = '{selected_lfc_gene}'
        AND pval_adj < 0.1
        ORDER BY ABS(lfc) DESC, pval_adj ASC
        LIMIT {num_regulators};
    """
    # Get LFCs for selected gene
    sig_regulators = query_database(query, CACHE_DIR, DATA_DIR)
    return sig_regulators.to_dict('records')

@app.callback(
    Output("table_sig_genes", "data"),
    Input("selected_lfc_target", "value"),
    Input("num_genes", "value"),
)  
def update_sig_genes_table(selected_lfc_target, num_genes):
    query = f"""
        SELECT Expressed_Gene_Symbol, lfc, pval_adj
        FROM lfc
        WHERE Target = '{selected_lfc_target}'
        AND pval_adj < 0.1
        ORDER BY ABS(lfc) DESC, pval_adj ASC
        LIMIT {num_genes};
    """
    # Get LFCs for selected target
    sig_genes = query_database(query, CACHE_DIR, DATA_DIR)
    return sig_genes.to_dict('records')

@app.callback(
    Output("target_cor_heatmap", "figure"),
    Output("target_cor_message", "children"),
    Input("selected_cor_targets", "value"),
)
def update_target_cor_heatmap(selected_cor_targets):
    if not selected_cor_targets:
        return go.Figure(), "Please select at least one target."
    if len(selected_cor_targets) > 25:
        return go.Figure(), "Only up to 25 targets can be visualized at once. Please select 1-25 targets for visualization."
    # Get correlations for selected targets
    filtered_target_cor = get_filtered_target_cor(selected_cor_targets, CACHE_DIR, DATA_DIR)
    if filtered_target_cor.empty:
        return go.Figure(), "No correlation data available for the selected targets."
    # Plot the correlation heatmap
    heatmap = create_target_cor_heatmap(filtered_target_cor, color_scale)
    return heatmap, ""

@app.callback(
    Output("target_scatter", "figure"),
    Input("selected_scatter_target1", "value"),
    Input("selected_scatter_target2", "value"),
)
def update_target_scatter(selected_scatter_target1, selected_scatter_target2):
    if not selected_scatter_target1 or not selected_scatter_target2:
        return go.Figure()
    query = f"""
        SELECT Target, Expressed_Gene_Symbol, lfc
        FROM lfc
        WHERE Target IN ('{selected_scatter_target1}', '{selected_scatter_target2}');
    """
    # Get LFCs for two selected target genes
    selected_lfc = query_database(query, CACHE_DIR, DATA_DIR)
    scatter_df = selected_lfc.pivot_table(index='Expressed_Gene_Symbol', columns='Target', values='lfc').reset_index()
    # Plot a scatterplot of the LFCs of one target vs. another one
    fig = px.scatter(scatter_df, x=selected_scatter_target1, y=selected_scatter_target2,
                     hover_name='Expressed_Gene_Symbol',
                     title=f"LFC after {selected_scatter_target1} vs {selected_scatter_target2} knockdown")
    fig.update_layout(template="plotly_white")
    fig.update_xaxes(title_text=f"LFC after {selected_scatter_target1} knockdown")
    fig.update_yaxes(title_text=f"LFC after {selected_scatter_target2} knockdown")
    return fig

@app.callback(
    Output("gene_cor_heatmap", "figure"),
    Output("gene_cor_message", "children"),
    Input("selected_cor_genes", "value"),
)
def update_gene_cor_heatmap(selected_cor_genes):
    if not selected_cor_genes:
        return go.Figure(), "Please select at least one gene."
    if len(selected_cor_genes) > 25:
        return go.Figure(), "Only up to 25 genes can be visualized at once. Please select 1-25 genes for visualization."
    # Get correlations for selected genes
    filtered_gene_cor = get_filtered_gene_cor(selected_cor_genes, CACHE_DIR, DATA_DIR)
    if filtered_gene_cor.empty:
        return go.Figure(), "No correlation data available for the selected genes."
    heatmap = create_gene_cor_heatmap(filtered_gene_cor, color_scale)
    return heatmap, ""

@app.callback(
    Output("gene_scatter", "figure"),
    Input("selected_scatter_gene1", "value"),
    Input("selected_scatter_gene2", "value"),
)
def update_gene_scatter(selected_scatter_gene1, selected_scatter_gene2):
    if not selected_scatter_gene1 or not selected_scatter_gene2:
        return go.Figure()
    query = f"""
        SELECT Target, Expressed_Gene_Symbol, lfc
        FROM lfc
        WHERE Expressed_Gene_Symbol IN ('{selected_scatter_gene1}', '{selected_scatter_gene2}');
    """
    # Get LFC for selected genes
    selected_lfc = query_database(query, CACHE_DIR, DATA_DIR)
    scatter_df = selected_lfc.pivot_table(index='Target', columns='Expressed_Gene_Symbol', values='lfc').reset_index()
    # Plot a scatterplot of the LFCs of one expressed gene vs. another one
    fig = px.scatter(scatter_df, x=selected_scatter_gene1, y=selected_scatter_gene2,
                     hover_name='Target',
                     title=f"LFC for {selected_scatter_gene1} vs {selected_scatter_gene2} expression")
    fig.update_layout(template="plotly_white")
    fig.update_xaxes(title_text=f"LFC for {selected_scatter_gene1} expression")
    fig.update_yaxes(title_text=f"LFC for {selected_scatter_gene2} expression")
    return fig

@app.callback(
    Output("cell_line_effects", "figure"),
    Input("selected_tar_target", "value"),
    Input("selected_tar_gene", "value"),
)    
def update_cell_line_effects(selected_tar_target, selected_tar_gene):
    selection = targeted_screen[
        (targeted_screen['Target'] == selected_tar_target) &
        (targeted_screen['Expressed_Gene_Symbol'] == selected_tar_gene)
    ].copy()
    if selection.empty:
        return go.Figure(layout=go.Layout(title="No data for selected target and gene"))
    selection['end'] = selection['wt_expr'] + selection['lfc']
    fig = create_heritability_plot(selection, selected_tar_target, selected_tar_gene)
    return fig

@app.callback(
    Output("table_heritable_gene", "data"),
    Input("selected_heritable_gene", "value"),
)
def update_heritable_gene_table(selected_heritable_gene):
    filtered_heritability = heritability[heritability['Downstream_Gene'] == selected_heritable_gene]
    return filtered_heritability.to_dict('records')

@app.callback(
    Output("table_heritable_target", "data"),
    Input("selected_heritable_target", "value"),
)
def update_heritable_target_table(selected_heritable_target):
    filtered_heritability = heritability[heritability['Target'] == selected_heritable_target]
    return filtered_heritability.to_dict('records')

if __name__ == '__main__':
    app.run(debug=True)

