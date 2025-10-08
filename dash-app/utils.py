import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import plotly.graph_objects as go
import json
import hashlib
import os

def get_filtered_lfc(targets, genes, CACHE_DIR, DATA_DIR):
    targets_sql = ', '.join(f"'{t}'" for t in targets)
    genes_sql = ', '.join(f"'{g}'" for g in genes)
    query = f"""
        SELECT Target, Expressed_Gene_Symbol, lfc 
        FROM lfc 
        WHERE Target IN ({targets_sql}) 
        AND Expressed_Gene_Symbol IN ({genes_sql})
    """
    df = query_database(query, CACHE_DIR, DATA_DIR)
    pivoted = df.pivot(index='Target', columns='Expressed_Gene_Symbol', values='lfc')
    return pivoted

def get_filtered_target_cor(targets, CACHE_DIR, DATA_DIR):
    targets_sql = ', '.join(f"'{t}'" for t in targets)
    query = f"""
        SELECT * 
        FROM target_cor 
        WHERE Target_1 IN ({targets_sql}) 
        AND Target_2 IN ({targets_sql})
    """
    df = query_database(query, CACHE_DIR, DATA_DIR)
    pivoted = df.pivot(index='Target_1', columns='Target_2', values='R')
    return pivoted

def get_filtered_gene_cor(genes, CACHE_DIR, DATA_DIR):
    genes_sql = ', '.join(f"'{g}'" for g in genes)
    query = f"""
        SELECT * 
        FROM gene_cor 
        WHERE Expressed_Gene_1 IN ({genes_sql}) 
        AND Expressed_Gene_2 IN ({genes_sql})
    """
    df = query_database(query, CACHE_DIR, DATA_DIR)
    pivoted = df.pivot(index='Expressed_Gene_1', columns='Expressed_Gene_2', values='R')
    return pivoted
    
def create_lfc_heatmap(data, colorscale='RdBu'):
    fig = go.Figure(data=go.Heatmap(
        z=data.values,
        x=data.columns,
        y=data.index,
        colorscale=colorscale,
        colorbar=dict(title='Log Fold Change'),
        zmid=0  # Center the color scale at 0
    ))
    fig.update_layout(
        title='Log Fold Change Heatmap',
        xaxis_title='Expressed Genes',
        yaxis_title='Targets',
    )
    return fig

def create_target_cor_heatmap(data, colorscale='RdBu'):
    fig = go.Figure(data=go.Heatmap(
        z=data.values,
        x=data.columns,
        y=data.index,
        colorscale=colorscale,
        colorbar=dict(title='Correlation'),
        zmin=-1, zmax=1
    ))
    fig.update_layout(
        title='Target Correlation Heatmap',
        xaxis_title='Target genes',
        yaxis_title='Target genes',
    )
    return fig

def create_gene_cor_heatmap(data, colorscale='RdBu'):
    fig = go.Figure(data=go.Heatmap(
        z=data.values,
        x=data.columns,
        y=data.index,
        colorscale=colorscale,
        colorbar=dict(title='Correlation'),
        zmin=-1, zmax=1
    ))
    fig.update_layout(
        title='Gene Correlation Heatmap',
        xaxis_title='Expressed Genes',
        yaxis_title='Expressed Genes',
    )
    return fig

def create_heritability_plot(selection, selected_tar_target, selected_tar_gene):
    
    fig = go.Figure()

    # Define a color map for individuals (still used for the plot lines and arrows)
    individual_colors = {ind: f'rgb({r},{g},{b})' for i, ind in enumerate(sorted(selection['individual'].unique())) for r, g, b in [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40), (148, 103, 189), (140, 86, 75), (227, 119, 194), (127, 127, 127), (188, 189, 34), (23, 190, 207)][i:i+1]}

    # Add segments and directional arrows for each individual
    for individual in sorted(selection['individual'].unique()):
        subset = selection[selection['individual'] == individual]
        color = individual_colors.get(individual, 'gray')
        for index, row in subset.iterrows():
            line_dash = 'solid' if row['linetype'] == 'solid' else 'dash'
            fig.add_trace(go.Scatter(
                x=[row['wt_expr'], row['end']],
                y=[row['Cell_Line'], row['Cell_Line']],
                mode='lines',
                line=dict(width=1.7, dash=line_dash, color=color),
                showlegend=False  # Don't show individual colors in the legend
            ))

            # Add directional arrowhead
            arrow_direction = 'arrow-right' if row['end'] > row['wt_expr'] else 'arrow-left'
            fig.add_trace(go.Scatter(
                x=[row['end']],
                y=[row['Cell_Line']],
                mode='markers',
                marker=dict(size=10, symbol=arrow_direction, color=color),
                showlegend=False  # Don't show individual colors in the legend
            ))

    # Add separate traces for the significance legend
    fig.add_trace(go.Scatter(
        x=[None],
        y=[None],
        mode='lines',
        line=dict(width=1.7, dash='solid', color='black'),
        name='Significant',
        legendgroup='significance',
        showlegend=True
    ))
    fig.add_trace(go.Scatter(
        x=[None],
        y=[None],
        mode='lines',
        line=dict(width=1.7, dash='dash', color='black'),
        name='Non-significant',
        legendgroup='significance',
        showlegend=True
    ))

    # Update layout
    fig.update_layout(
        title=f"Effects of {selected_tar_target} knockdown on {selected_tar_gene} expression",
        yaxis_title="Cell Line",
        xaxis_title="Log expression before and after knockdown",
        template="plotly_white",
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="right",
            x=1,
            tracegroupgap=20,
            traceorder='normal' # Ensure significance legend appears last
        )
    )
    return fig

def generate_cache_key(query):
    """Generates an MD5 hash based on the query and its parameters."""
    combined = query.encode('utf-8')
    return hashlib.md5(combined).hexdigest() + '.json'

def get_cached_data(CACHE_DIR, cache_key):
    """Retrieves data from the cache if the file exists and is valid."""
    cache_path = os.path.join(CACHE_DIR, cache_key)
    if os.path.exists(cache_path):
        try:
            with open(cache_path, 'r') as f:
                return json.load(f)
        except json.JSONDecodeError:
            os.remove(cache_path)  # Remove corrupted cache file
            return None
    return None

def save_data_to_cache(data, CACHE_DIR, cache_key):
    """Saves data to the cache as a JSON file."""
    cache_path = os.path.join(CACHE_DIR, cache_key)
    with open(cache_path, 'w') as f:
        json.dump(data, f)

def query_database(query, CACHE_DIR, DATA_DIR):
    """Queries the SQLite database and returns the results, with caching."""
    cache_key = generate_cache_key(query)
    cached_data = get_cached_data(CACHE_DIR, cache_key)

    if cached_data:
        return pd.DataFrame(cached_data)  

    conn = sqlite3.connect(DATA_DIR + 'data.sqlite')
    df = pd.read_sql_query(query, conn)
    conn.close()
    data_to_cache = df.to_dict('records')
    save_data_to_cache(data_to_cache, CACHE_DIR, cache_key)
    return df