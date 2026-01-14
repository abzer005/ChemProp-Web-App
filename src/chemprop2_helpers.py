import streamlit as st
import os
import zipfile
import tempfile
from io import BytesIO
from .common import *  # Importing common functionalities from the 'common' module
from src.cleanup import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.stats import spearmanr, pearsonr
import plotly.graph_objects as go
import networkx as nx
from streamlit_agraph import agraph, Node, Edge, Config
from src.chemprop1_helpers import _legend_html, render_spectra_modifinder
from src.chemprop1_helpers import drop_blank_score_rows
from src.chemprop1_helpers import add_mz_rt_from_ft
from src.chemprop1_helpers import add_names_from_gnps
from src.chemprop1_helpers import render_spectra_modifinder


def get_chemprop2_inputs_from_state():
    required = ["nw", "chemprop_ft", "chemprop_md"]
    if not all(k in st.session_state for k in required):
        return None

    network_df = st.session_state["nw"].copy()
    features_df = st.session_state["chemprop_ft"].copy()
    metadata_df = st.session_state["chemprop_md"].copy()

    if network_df.empty or features_df.empty or metadata_df.empty:
        return None

    return network_df, features_df, metadata_df

def chemprop2_controls():
    if "run_chemprop2" not in st.session_state:
        st.session_state.run_chemprop2 = False

    show_options = st.checkbox("Run ChemProp2 score calculation", value=False)
    st.session_state.run_chemprop2 = bool(show_options)

    if not show_options:
        return False, None

    mode = st.radio(
        "Run ChemProp2 on:",
        ("Provided Edge Table", "Cascade edges", "User defined edge"),
        horizontal=True,
        key="chemprop2_mode"
    )

    return True, mode


# drop_blank_score_rows, add_mz_rt_from_ft, add_names_from_gnps are imported from chemprop1_helpers.py

def run_chemprop2_pipeline(network_df, features_df, metadata_df):
    
    chemprop2_df = run_chemProp2_scoring(network_df, features_df, metadata_df)

    if chemprop2_df is None or chemprop2_df.empty:
        return None

    chemprop2_df = drop_blank_score_rows(chemprop2_df, base_cols=network_df.shape[1])
    chemprop2_df = add_mz_rt_from_ft(chemprop2_df, st.session_state["ft"])
    
    # Add names if GNPS annotations exist
    if "an_gnps" in st.session_state and isinstance(st.session_state["an_gnps"], pd.DataFrame):
        gnps_df = st.session_state["an_gnps"]
        if not gnps_df.empty:
            chemprop2_df = add_names_from_gnps(chemprop2_df, gnps_df)

    st.session_state["ChemProp2_scores"] = chemprop2_df
    st.write(
        f"ChemProp2 Scoring Results " 
        f"({chemprop2_df.shape[0]} rows × {chemprop2_df.shape[1]} columns):"
        )
    st.dataframe(chemprop2_df, hide_index=True, use_container_width=True)
    return chemprop2_df

def run_chemProp2_scoring(desired_network_table, desired_feature_table, chemprop_metadata):

    names_nodes = list(desired_feature_table.index)
    
    # Initialize the score columns with NA values (using np.nan for missing values)
    ChemProp2 = np.full(len(desired_network_table), np.nan)
    ChemProp_spearman = np.full(len(desired_network_table), np.nan)
    ChemProp_log = np.full(len(desired_network_table), np.nan)
    ChemProp_sqrt = np.full(len(desired_network_table), np.nan)
    
    for i in range(len(desired_network_table)):
        try:
            # Attempt to find the indices of the feature IDs
            feature_ids = [names_nodes.index(desired_network_table.iloc[i, 0]), 
                           names_nodes.index(desired_network_table.iloc[i, 1])]
            
            # Subset feature data
            x = desired_feature_table.iloc[feature_ids, :]
            
            # Get metadata timepoints
            chemprop_sample_names = chemprop_metadata.index
            chemprop_timepoint = chemprop_metadata.iloc[:, 0]
            
            # Reorder x based on sample names in metadata
            reorder_id = [list(desired_feature_table.columns).index(name) for name in chemprop_sample_names]
            reordered_x = x.iloc[:, reorder_id].T
            reordered_x = pd.concat([chemprop_timepoint.reset_index(drop=True), reordered_x.reset_index(drop=True)], axis=1)
            
            # Pearson correlation
            corr_result = reordered_x.corr(method="pearson")
            ChemProp2[i] = (corr_result.iloc[0, 2] - corr_result.iloc[0, 1]) / 2
            
            # Spearman correlation
            corr_spearman = reordered_x.corr(method="spearman")
            ChemProp_spearman[i] = (corr_spearman.iloc[0, 2] - corr_spearman.iloc[0, 1]) / 2
            
            # Log transformation
            log_reorderedX = pd.concat([reordered_x.iloc[:, 0], np.log1p(reordered_x.iloc[:, 1:3])], axis=1)
            corr_log = log_reorderedX.corr(method="pearson")
            ChemProp_log[i] = (corr_log.iloc[0, 2] - corr_log.iloc[0, 1]) / 2
            
            # Square root transformation
            sqrt_reorderedX = pd.concat([reordered_x.iloc[:, 0], np.sqrt(reordered_x.iloc[:, 1:3])], axis=1)
            corr_sqrt = sqrt_reorderedX.corr(method="pearson")
            ChemProp_sqrt[i] = (corr_sqrt.iloc[0, 2] - corr_sqrt.iloc[0, 1]) / 2
        
        except ValueError:
            # If the feature ID is not found, leave the ChemProp scores as np.nan (NA in Pandas)
            ChemProp2[i] = np.nan
            ChemProp_spearman[i] = np.nan
            ChemProp_log[i] = np.nan
            ChemProp_sqrt[i] = np.nan

    # Combine results into a new DataFrame
    new_network_table = desired_network_table.copy()
    new_network_table['ChemProp2'] = ChemProp2
    new_network_table['ChemProp_spearman'] = ChemProp_spearman
    new_network_table['ChemProp_log'] = ChemProp_log
    new_network_table['ChemProp_sqrt'] = ChemProp_sqrt

    # Add absolute value columns
    abs_values = new_network_table[['ChemProp2', 'ChemProp_spearman', 'ChemProp_log', 'ChemProp_sqrt']].abs()
    abs_values.columns = [f"abs_{col}" for col in abs_values.columns]
    abs_values = abs_values.apply(pd.to_numeric, errors='coerce')
    
    # Add sign of ChemProp2
    new_network_table['Sign_ChemProp2'] = np.sign(new_network_table['ChemProp2'])
    
    # Concatenate absolute values
    ChemProp2_file = pd.concat([new_network_table, abs_values], axis=1)

    return ChemProp2_file


def render_chemprop2_results_summary(network_df, result_df):

    st.info(
        f"""
    **Edge Filtering Summary**

    - 📊 **Original edge table:** {network_df.shape[0]} edges  
    - 🧹 **ChemProp2 edge table:** {result_df.shape[0]} edges
    """
    )

    if network_df.shape[0] != result_df.shape[0]:
        st.warning("The reduced number of edges might be due to the removal of blank entries.")


def render_chemprop2_download_buttons(result_df):
    user_filename = st.text_input(
        "Enter the filename for the CSV and press Enter to apply the name:",
        value="chemprop2_scores_results.csv",
    )
    if not user_filename.endswith(".csv"):
        user_filename += ".csv"

    data_csv = result_df.to_csv(index=False).encode("utf-8")

    col1, col2, col3, _ = st.columns([1.1, 1, 1, 4.9])

    with col1:
        st.download_button(
            label="Download Results as CSV",
            data=data_csv,
            file_name=user_filename,
            mime="text/csv",
        )

    with col2:
        # remove the __name__ guard in Streamlit; it’s not needed
        if st.button("Download GraphML (ZIP)"):
            generate_graphml_zip_chemprop2()

    with col3:
        gnps_task_id = st.session_state.get("gnps_task_id")

        if gnps_task_id:
            gnps_url = (
                "https://www.gnps2.org/dashboards/networkviewer/"
                f"?usi=mzdata%3AGNPS2%3ATASK-{gnps_task_id}-nf_output%2Fnetworking%2Fnetwork.graphml"
                f"&usi-mgf=mzdata%3AGNPS2%3ATASK-{gnps_task_id}-nf_output%2Fclustering%2Fspectra_reformatted.mgf#{{}}"
            )

            st.link_button(
                "Visualize Network in GNPS2",
                gnps_url,
                type="primary",
            )
        else:
            st.info("ℹ️ GNPS task ID not available for network visualization.")



def generate_graphml_zip_chemprop2():
    if 'ChemProp2_scores' in st.session_state:
        df = st.session_state['ChemProp2_scores'].copy()

        chemprop2_graphml_filename = "chemprop_graph.graphml"
        style_filename = "resources/ChemProp2_styles.xml"

        with tempfile.TemporaryDirectory() as temp_dir:
            chemprop2_graphml_path = os.path.join(temp_dir, chemprop2_graphml_filename)
            generate_graphml_with_secondary_edges_chemprop2(df, chemprop2_graphml_path)

            # Full path to styles.xml
            style_path = os.path.join(os.path.dirname(__file__), style_filename)
            
            # Check if the styles.xml file exists at the calculated path
            if not os.path.exists(style_path):
                st.error("The styles.xml file is not found in the 'resources' folder.")
                return

            # Create the zip file in the temporary directory
            zip_path = os.path.join(temp_dir, 'chemprop2_graph_with_styles.zip')
            with zipfile.ZipFile(zip_path, 'w') as zipf:
                zipf.write(chemprop2_graphml_path, arcname=chemprop2_graphml_filename)
                zipf.write(style_path, arcname="styles.xml")

            # Provide the zip file for download
            with open(zip_path, 'rb') as f:
                st.download_button(
                    label="Download GraphML and Styles as ZIP",
                    data=f,
                    file_name='chemprop2_graph_with_styles.zip',
                    mime='application/zip'
                )

def render_chemprop2_scores_plot(scores_df):
    df = scores_df.copy()
    df["CLUSTER IDs"] = df["CLUSTERID1"].astype(str) + "-" + df["CLUSTERID2"].astype(str)

    fig = px.scatter(
        df,
        x="DeltaMZ",
        y="ChemProp2",
        hover_name="CLUSTER IDs",
        title="Scatter Plot: DeltaMZ vs ChemProp2 scores",
        labels={"DeltaMZ": "Delta M/Z", "ChemProp2": "ChemProp2 scores"},
    )
    st.plotly_chart(fig, use_container_width=True)


def plot_intensity_trends_single_row_chemprop2(row, feature_table, metadata):
    """
    ChemProp2 intensity trend (scatter + mean dashed lines) for a single edge row.

    Returns:
        fig (plotly.graph_objects.Figure)
        status (dict): info about missing data / issues
    """
    status = {
        "missing_ids": [],
        "no_overlap_samples": False,
        "metadata_missing": False,
    }

    clusterid1 = row.get("CLUSTERID1", None)
    clusterid2 = row.get("CLUSTERID2", None)

    fig = go.Figure()
    fig.update_layout(template="plotly_white", width=800, height=600)

    if clusterid1 is None or clusterid2 is None:
        status["missing_ids"] = ["CLUSTERID1/CLUSTERID2"]
        fig.add_annotation(text="Missing CLUSTER IDs.", showarrow=False)
        return fig, status

    # ---- normalize feature table index and cluster ids safely
    ft = feature_table.copy()
    ft.index = ft.index.astype(str)

    def _to_id_str(x):
        """Robustly convert 123 / '123' / 123.0 / '123.0' -> '123'."""
        try:
            return str(int(float(x)))
        except Exception:
            return str(x)

    id1 = _to_id_str(clusterid1)
    id2 = _to_id_str(clusterid2)

    if id1 not in ft.index:
        status["missing_ids"].append(clusterid1)
    if id2 not in ft.index:
        status["missing_ids"].append(clusterid2)

    if status["missing_ids"]:
        fig.add_annotation(
            text=f"Feature(s) not found in feature table: {status['missing_ids']}",
            xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False
        )
        return fig, status

    if metadata is None or metadata.empty:
        status["metadata_missing"] = True
        fig.add_annotation(text="Metadata is missing.", showarrow=False)
        return fig, status

    time_col = "Timepoint" if "Timepoint" in metadata.columns else metadata.columns[0]

    common_samples = ft.columns.intersection(metadata.index)
    if len(common_samples) == 0:
        status["no_overlap_samples"] = True
        fig.add_annotation(
            text="No overlapping samples between feature table and metadata.",
            showarrow=False
        )
        return fig, status

    # ---- build plotting dataframe (numeric coercion prevents weird dtype issues)

    label1 = f"ID {clusterid1}"
    label2 = f"ID {clusterid2}"

    df = pd.DataFrame({
        "Timepoint": metadata.loc[common_samples, time_col].astype(str).values,
        label1: pd.to_numeric(ft.loc[id1, common_samples], errors="coerce").values,
        label2: pd.to_numeric(ft.loc[id2, common_samples], errors="coerce").values,
    })

    # long format for points
    long_df = df.melt(
        id_vars="Timepoint",
        value_vars=[label1, label2],
        var_name="Feature",
        value_name="Value"
    )

    # ---- fixed color mapping (scatter + means)
    color_map = {
        label1: "blue",  # blue
        label2: "#d62728",  # red
    }

    # ---- scatter points (per sample) + mean dashed lines (per timepoint)
    fig = px.scatter(
        long_df,
        x="Timepoint",
        y="Value",
        color="Feature",
        color_discrete_map=color_map,
        template="plotly_white",
        title=f"Intensity Trend Scatter Plot {clusterid1} vs {clusterid2}",
        labels={"Timepoint": "Time", "Value": "Feature Intensity"},
        width=800,
        height=600,
    )

    means = (
        df.groupby("Timepoint")[[label1, label2]]
        .mean(numeric_only=True)
        .reset_index()
    )

    fig.add_scatter(
        x=means["Timepoint"],
        y=means[label1],
        mode="lines",
        name=f"Mean {label1}",
        line=dict(color=color_map[label1], dash="dash"),
    )
    fig.add_scatter(
        x=means["Timepoint"],
        y=means[label2],
        mode="lines",
        name=f"Mean {label2}",
        line=dict(color=color_map[label2], dash="dash"),
    )

    # ---- ChemProp2 annotation (score + names)
    score = row.get("ChemProp2", np.nan)
    name1 = row.get("ID1_name", "")
    name2 = row.get("ID2_name", "")

    fig.add_annotation(
        text=f"ChemProp2 score = {'NA' if pd.isna(score) else f'{score:.2f}'}<br>"
             f"{name1} <b>VS</b> {name2}",
        xref="paper", yref="paper",
        x=0.5, y=1.1,
        showarrow=False,
    )

    return fig, status


def generate_graphml_with_secondary_edges_chemprop2(df, output_file):
    """
    Generate a GraphML file from the ChemProp2_scores DataFrame
    
    Parameters:
    df (pd.DataFrame): The ChemProp2_scores DataFrame containing node and edge data.
    output_file (str): The name of the GraphML file to generate.
    """
    # Create a directed graph allowing multiple edges
    G = nx.MultiDiGraph()

    # Check for 'chemprop_ft' in session state for additional node data
    node_table = st.session_state['chemprop_ft'] if 'chemprop_ft' in st.session_state else None

    # Add nodes and edges with attributes
    for _, row in df.iterrows():
        clusterid1 = row['CLUSTERID1']
        clusterid2 = row['CLUSTERID2']

        # Extract node names
        id1_name = row['ID1_name'] if 'ID1_name' in row else str(clusterid1)
        id2_name = row['ID2_name'] if 'ID2_name' in row else str(clusterid2)
        
        # Add nodes with optional attributes
        G.add_node(clusterid1, node_names=id1_name)
        G.add_node(clusterid2, node_names=id2_name)

        # Add additional attributes from node_table if available
        if node_table is not None:
            if clusterid1 in node_table.index:
                for col in node_table.columns:
                    G.nodes[clusterid1][col] = node_table.loc[clusterid1, col]
            if clusterid2 in node_table.index:
                for col in node_table.columns:
                    G.nodes[clusterid2][col] = node_table.loc[clusterid2, col]

        # Define edge attributes for the primary edge
        primary_edge_attributes = {
            'ComponentIndex': row['ComponentIndex'],
            'Cosine': row['Cosine'],
            'DeltaMZ': row['DeltaMZ'],
            'ChemProp2': row['ChemProp2'],
            'ChemProp_spearman': row['ChemProp_spearman'],
            'ChemProp_log': row['ChemProp_log'],
            'ChemProp_sqrt': row['ChemProp_sqrt'],
            'label': f"∆mz: {row['DeltaMZ']:.2f}",
            'color': "black"
        }

        # Add primary edge
        G.add_edge(clusterid1, clusterid2, key="primary", **primary_edge_attributes)

        # Define secondary edge attributes based on abs_ChemProp2 and Sign_ChemProp2
        abs_chemprop = row['abs_ChemProp2']
        sign_chemprop2 = row['Sign_ChemProp2']
        secondary_edge_attributes = {
                'weight': abs_chemprop,
                'label': f"{abs_chemprop:.2f}",
                'color': "red",
                'abs_ChemProp2': row['abs_ChemProp2'],
                'abs_ChemProp_spearman': row['abs_ChemProp_spearman'],
                'abs_ChemProp_log': row['abs_ChemProp_log'],
                'abs_ChemProp_sqrt': row['abs_ChemProp_sqrt']
                }
        # Use Sign_ChemProp2 to determine direction of the secondary edge
        if sign_chemprop2 == 1:
            G.add_edge(clusterid1, clusterid2, key="secondary", **secondary_edge_attributes, arrow=True)
        elif sign_chemprop2 == -1:
            G.add_edge(clusterid2, clusterid1, key="secondary", **secondary_edge_attributes, arrow=True)

    # Write the graph to GraphML
    nx.write_graphml(G, output_file)

    return output_file


def calculate_fdr(edge_df, features_df, metadata_df, score_range=(-1, 1), bin_size=0.001):
    """
    Calculates a cutoff to set a false discovery rate (FDR) for correlation scores. A target-decoy approach 
    is used to calculate FDR scores, allowing users to select a range to trust the ChemProp2 scores.

    Parameters:
    - edge_df (pd.DataFrame): The DataFrame containing node and edge data.
    - features_df (pd.DataFrame): Feature table.
    - metadata_df (pd.DataFrame): Metadata table.
    - score_range (tuple): Range for score binning (default is (-1, 1)).
    - bin_size (float): Size of bins for histogram (default is 0.01).
    Returns:
    - combined_fdr_df (pd.DataFrame): DataFrame with FDR values across bins for the decoy set.
    """
    # Calculate target scores
    target_df = run_chemProp2_scoring(edge_df, features_df, metadata_df)

    # Define bins based on specified range and bin size
    bins = np.arange(score_range[0], score_range[1] + bin_size, bin_size)

    # get the decoy table
    decoy_features_df = features_df.apply(np.random.permutation)
    decoy_features_df += np.random.normal(0, 0.1, decoy_features_df.shape)  # Add noise
        
    # decoy_edge_df = edge_df.copy()
    # decoy_edge_df['CLUSTERID1'] = np.random.permutation(edge_df['CLUSTERID1'])
    # decoy_edge_df['CLUSTERID2'] = np.random.permutation(edge_df['CLUSTERID2'])

    decoy_df = run_chemProp2_scoring(edge_df, decoy_features_df, metadata_df)

    # Calculate target and decoy counts per bin
    target_counts, _ = np.histogram(target_df['ChemProp2'], bins=bins)
    decoy_counts, _ = np.histogram(decoy_df['ChemProp2'], bins=bins)
    
    # Calculate FDR for each bin
    fdr_df = pd.DataFrame({'Range_min': bins[:-1], 
                               'Range_max': bins[1:], 
                               'Target_counts': target_counts,
                               'Decoy_counts': decoy_counts
                               })

    epsilon = 0

    # Separate positive and negative bins
    positive_bins = fdr_df[fdr_df['Range_max'] > 0].reset_index(drop=True)
    negative_bins = fdr_df[fdr_df['Range_max'] < 0].reset_index(drop=True)

    # Calculate cumulative FDR for positive bins (0 to 1, forward order)
    positive_bins = positive_bins[::-1].reset_index(drop=True)  # Reverse order for cumulative sum
    positive_bins['FDR'] = positive_bins['Decoy_counts'].cumsum() / (
        positive_bins['Target_counts'].cumsum()  + epsilon + positive_bins['Decoy_counts'].cumsum()
    )

    positive_bins = positive_bins[::-1].reset_index(drop=True)  # Restore original order

    # Calculate cumulative FDR for negative bins (-1 to 0)
    negative_bins['FDR'] = negative_bins['Decoy_counts'].cumsum() / (
        negative_bins['Target_counts'].cumsum() +  epsilon + negative_bins['Decoy_counts'].cumsum()
    )
   
    # Combine positive and negative bins
    combined_fdr_df = pd.concat([negative_bins, positive_bins]).reset_index(drop=True)

    # Plot FDR results
    fig_fdr = go.Figure()
    fig_fdr.add_trace(go.Scatter(
        x=combined_fdr_df['Range_min'],
        y=combined_fdr_df['FDR']*100,
        mode='markers',
        name = 'FDR Score'
    ))

    # Define FDR threshold levels and their colors
    fdr_thresholds = [1, 5, 10]
    colors = ['green', 'orange', 'red']

    ########################################
    # Sort the combined_fdr_df by FDR values
    combined_fdr_df_sorted = combined_fdr_df.sort_values(by='FDR').reset_index(drop=True)
    combined_fdr_df_sorted['FDR'] *= 100

    # Loop through each threshold to add a horizontal line and highlight the closest point
    for fdr_value, color in zip(fdr_thresholds, colors):
        # Add a horizontal line for the FDR threshold
        fig_fdr.add_shape(
            type="line",
            x0=combined_fdr_df['Range_min'].min(), 
            x1=combined_fdr_df['Range_min'].max(),
            y0=fdr_value, y1=fdr_value,
            line=dict(color=color, width=1, dash="dash"),
            name=f"{fdr_value}% FDR"
        )

        # Identify points that are within a small tolerance of the FDR threshold
        tolerance = 0.8  # Adjust tolerance as needed to capture points near each threshold
        threshold_points = combined_fdr_df[(combined_fdr_df['FDR'] * 100 >= fdr_value - tolerance) & 
                                        (combined_fdr_df['FDR'] * 100 <= fdr_value + tolerance)]

        # Add a marker at the closest point to the FDR threshold
        fig_fdr.add_trace(go.Scatter(
            x=threshold_points['Range_min'],
            y=threshold_points['FDR'] * 100,
            mode='markers',
            marker=dict(color=color, size=10, symbol='circle-open'),
            name=f"{fdr_value}% Threshold"
        ))

    # Update layout for the figure
    fig_fdr.update_layout(
        title="Overlay of FDR Score for Target-Decoy Sets",
        xaxis_title="ChemProp2 Score Range",
        yaxis_title="FDR (%)",
        width=900,
        height=500
    )

    # Plot histograms for target and decoys
    fig_histogram = go.Figure()
    fig_histogram.add_trace(go.Histogram(
        x=target_df['ChemProp2'],
        name='Target',
        xbins=dict(start=-1, end=1, size=0.1),
        opacity=0.5,
        marker_color='blue'
    ))

    fig_histogram.add_trace(go.Histogram(
        x=decoy_df['ChemProp2'],
        name='Decoy',
        xbins=dict(start=-1, end=1, size=0.1),
        opacity=0.3,
        marker_color='red'
    ))

    fig_histogram.update_layout(
        title="Histogram of Target vs. Decoy Scores",
        xaxis_title="ChemProp2 Score Range",
        yaxis_title="Frequency",
        xaxis_range=[-1, 1],
        barmode='overlay',
        width=900,
        height=500
    )
    
    return combined_fdr_df, fig_histogram, fig_fdr


def render_chemprop2_filters_and_plots(scores_df, features_df, metadata_df):
    if scores_df is None or scores_df.empty:
        return

    if not st.checkbox("Show Filters", key="show_filters_chemprop2"):
        return
    
    if scores_df.shape[0] == 1:
        filtered_df = scores_df.copy()
    
    else:
        df = scores_df.copy()
        df["CLUSTER IDs"] = df["CLUSTERID1"].astype(str) + "-" + df["CLUSTERID2"].astype(str)

        filter_mode = st.radio(
            "Filter type:",
            ("Filter by Score", "Filter by M/Z Range", "Filter by Name", "Filter by Cluster ID"),
            horizontal=True,
        )

        filtered_df = get_chemprop2_filtered_edges_ui(df, filter_mode)
        st.dataframe(filtered_df, hide_index=True, use_container_width=True)
    render_chemprop2_edge_detail_plots(filtered_df, scores_df, features_df, metadata_df)
    if (st.session_state.get("gnps_task_id") or "").strip():
        render_spectra_modifinder(filtered_df, st.session_state.get("an_gnps"))


def get_chemprop2_filtered_edges_ui(df, filter_mode):
    c1, c2, c3 = st.columns(3)

    if filter_mode == "Filter by Score":
        with c1:
            score_min = st.number_input("Min Score:", value=float(df["ChemProp2"].min()))
        with c2:
            score_max = st.number_input("Max Score:", value=float(df["ChemProp2"].max()))

        sub = df[(df["ChemProp2"] >= score_min) & (df["ChemProp2"] <= score_max)].copy()
        with c3:
            if sub.empty:
                st.warning("No edges in this score range.")
                return sub
            sub["Dropdown_Display"] = (
                "ID: " + sub["CLUSTER IDs"].astype(str) + ", Score: " + sub["ChemProp2"].round(3).astype(str)
            )
            choice = st.selectbox("Select the edge to view plots", sub["Dropdown_Display"].tolist())
            return sub[sub["Dropdown_Display"] == choice].drop(columns=["Dropdown_Display"])

    if filter_mode == "Filter by M/Z Range":
        with c1:
            mz_min = st.number_input("Min DeltaMZ:", value=float(df["DeltaMZ"].min()))
        with c2:
            mz_max = st.number_input("Max DeltaMZ:", value=float(df["DeltaMZ"].max()))

        sub = df[(df["DeltaMZ"] >= mz_min) & (df["DeltaMZ"] <= mz_max)].copy()

        with c3:
            if sub.empty:
                st.warning("No edges in this DeltaMZ range.")
                return sub
            sub["Dropdown_Display"] = (
                "ID: " + sub["CLUSTER IDs"].astype(str) + ", DeltaMZ: " + sub["DeltaMZ"].round(3).astype(str)
            )
            choice = st.selectbox("Select the edge to view plots", sub["Dropdown_Display"].tolist())
            return sub[sub["Dropdown_Display"] == choice].drop(columns=["Dropdown_Display"])

    if filter_mode == "Filter by Name":
        if not {"ID1_name", "ID2_name"}.issubset(df.columns):
            st.warning("You don't have names for the cluster IDs in the table.")
            return df.iloc[0:0]  # empty

        with c1:
            text = st.text_input("Enter text:", "")
        sub = df.copy()
        if text:
            sub = sub[
                sub["ID1_name"].str.contains(text, case=False, na=False)
                | sub["ID2_name"].str.contains(text, case=False, na=False)
            ].copy()

        with c2:
            if sub.empty:
                st.warning("No edges match that text.")
                return sub
            choice = st.selectbox("Select the edge to view plots", sub["CLUSTER IDs"].tolist())
            return sub[sub["CLUSTER IDs"] == choice]

    # Filter by Cluster ID
    choice = st.selectbox("Select the edge to view plots:", options=df["CLUSTER IDs"].unique())
    return df[df["CLUSTER IDs"] == choice]


def render_chemprop2_edge_detail_plots(filtered_df, all_scores_df, features_df, metadata_df):
    if filtered_df is None or filtered_df.empty:
        return

    with st.container():
        c1, c2 = st.columns(2)

        with c1:
            selected_row = filtered_df.iloc[0]
            fig, status = plot_intensity_trends_single_row_chemprop2(
                selected_row,
                features_df,
                metadata_df
            )

            # ---- Streamlit warnings ----
            if status["missing_ids"]:
                st.warning(
                    f"Cannot plot intensities for feature(s) not present in the feature table: "
                    f"{status['missing_ids']}"
                )

            if status["metadata_missing"]:
                st.warning("Metadata table is missing or empty. Intensity trends cannot be plotted.")

            if status["no_overlap_samples"]:
                st.warning(
                    "No overlapping samples between feature table and metadata. "
                    "Please check sample naming consistency."
                )

            st.plotly_chart(fig, use_container_width=True)

        with c2:
            try:
                comp_index = filtered_df["ComponentIndex"].iloc[0]
                plot_df = all_scores_df[all_scores_df["ComponentIndex"] == comp_index].copy()
                allowed = ['ComponentIndex', 'Cosine', 'DeltaMZ', 'ChemProp2', 
                           'ChemProp_spearman', 'ChemProp_log', 'ChemProp_sqrt', 
                           'Sign_ChemProp2', 'abs_ChemProp2', 'abs_ChemProp_spearman', 
                           'abs_ChemProp_log', 'abs_ChemProp_sqrt'
                           ]
                available = [c for c in plot_df.columns if c in allowed]

                edge_label_column = st.selectbox(
                    "Select column for edge labels:",
                    options=available,
                    help="To save the network image as PNG, right-click on empty space and select 'Save image as'.",
                )

                nodes, edges = generate_graph_from_df_chemprop2(plot_df, filtered_df, edge_label_column)

                config = Config(
                    width=800, height=600,
                    directed=True,
                    nodeHighlightBehavior=True,
                    highlightColor="#F7A7A6",
                    collapsible=True,
                    node={"labelProperty": "label"},
                    link={"labelProperty": "label", "renderLabel": True},
                    staticGraph=False,
                )

                st.markdown(_legend_html(), unsafe_allow_html=True)
                agraph(nodes=nodes, edges=edges, config=config)

            except ModuleNotFoundError:
                st.error("This page requires the `pygraphviz` package, which is not available in the Windows app.")

def generate_graph_from_df_chemprop2(df: pd.DataFrame, filtered_df: pd.DataFrame, edge_label_column: str):
    """
    Build nodes/edges for streamlit-agraph from ChemProp2 score dataframe.

    - Creates one node per unique cluster id in df (and ensures source/target are present)
    - Node label: m/z (2 decimals) when available, otherwise the cluster id
    - Node title: id, name, m/z, RT (best-effort)
    - Colors: source=blue, target=red, others=lightgray
    - Edges: direction based on Sign_ChemProp2 when abs_ChemProp2 > 0
    """

    nodes = []
    edges = []
    added = set()

    if filtered_df is None or filtered_df.empty:
        return nodes, edges

    # Normalize source/target ids as strings (agraph node ids are strings)
    source_id = str(filtered_df["CLUSTERID1"].iloc[0])
    target_id = str(filtered_df["CLUSTERID2"].iloc[0])

    # Work on a copy; normalize id columns to string for consistent comparisons
    df2 = df.copy()
    df2["CLUSTERID1"] = df2["CLUSTERID1"].astype(str)
    df2["CLUSTERID2"] = df2["CLUSTERID2"].astype(str)

    # Build neighbor map (cluster id -> neighbor cluster id)
    neighbor_map = {}
    if "Neighbor" in df2.columns:
        tmp = df2[["CLUSTERID2", "Neighbor"]].dropna()
        tmp["CLUSTERID2"] = tmp["CLUSTERID2"].astype(str)

        # If duplicates exist, keep the smallest neighbor distance
        tmp["Neighbor"] = pd.to_numeric(tmp["Neighbor"], errors="coerce")
        neighbor_map = tmp.groupby("CLUSTERID2")["Neighbor"].min().to_dict()

    # Helper: add node once with best-available metadata
    def add_node(node_id: str, mz=None, rt=None, name=None, color="lightgray", neighbor=None):
        if node_id in added:
            return

        label = f"{mz:.2f}" if pd.notna(mz) else node_id

        title_parts = [f"ID: {node_id}"]
        if pd.notna(name) and str(name).strip().lower() != "nan":
            title_parts.append(f"Name: {name}")
        if pd.notna(mz):
            title_parts.append(f"m/z: {mz}")
        if pd.notna(rt):
            title_parts.append(f"RT: {rt}")

        # only add if neighbor was passed
        if neighbor is not None and pd.notna(neighbor):
            title_parts.append(f"Neighbor: {int(neighbor)}")
        title = "\n".join(title_parts)

        nodes.append(Node(
            id=node_id,
            label=label,
            size=20,
            color=color,
            title=title
        ))
        added.add(node_id)

    # 1) Build a metadata lookup for any node id from either side of edges
    #    First occurrence wins (good enough for mz/rt/name)
    meta = {}
    for _, r in df2.iterrows():
        a = r["CLUSTERID1"]
        b = r["CLUSTERID2"]

        # Side A metadata
        if a not in meta:
            meta[a] = {
                "mz": r.get("ID1_mz", pd.NA),
                "rt": r.get("ID1_RT", pd.NA),
                "name": r.get("ID1_name", pd.NA),
            }
        # Side B metadata
        if b not in meta:
            meta[b] = {
                "mz": r.get("ID2_mz", pd.NA),
                "rt": r.get("ID2_RT", pd.NA),
                "name": r.get("ID2_name", pd.NA),
            }

    # Ensure source/target are present in meta (even if not found)
    meta.setdefault(source_id, {"mz": pd.NA, "rt": pd.NA, "name": pd.NA})
    meta.setdefault(target_id, {"mz": pd.NA, "rt": pd.NA, "name": pd.NA})

    # 2) Add nodes (all unique ids)
    # Color rule: source blue, target red, others lightgray
    for node_id, m in meta.items():
        color = "blue" if node_id == source_id else ("red" if node_id == target_id else "lightgray")
        nb = neighbor_map.get(node_id, None)
        add_node(node_id, mz=m.get("mz"), rt=m.get("rt"), name=m.get("name"), color=color, neighbor=nb)

    # 3) Add edges
    #    - label = selected edge_label_column (best-effort formatting)
    #    - direction based on Sign_ChemProp2 (1: a->b, -1: b->a)
    for _, r in df2.iterrows():
        a = r["CLUSTERID1"]
        b = r["CLUSTERID2"]

        abs_cp = r.get("abs_ChemProp2", 0)
        sign = r.get("Sign_ChemProp2", 0)

        # Skip if abs score is missing/zero/NaN
        if pd.isna(abs_cp) or float(abs_cp) <= 0:
            continue

        # Edge label
        val = r.get(edge_label_column, "")
        if pd.isna(val):
            label = ""
        else:
            # try numeric formatting, otherwise string
            try:
                label = f"{float(val):.2f}"
            except Exception:
                label = str(val)

        if sign == 1:
            src, tgt = a, b
        elif sign == -1:
            src, tgt = b, a
        else:
            # if sign missing, default a->b
            src, tgt = a, b

        edges.append(Edge(
            source=src,
            target=tgt,
            label=label,
            color="orange",
            arrow=True,
            font={"size": 10},
        ))

    return nodes, edges

# def generate_graph_from_df_chemprop2(df, filtered_df, edge_label_column):
#     """
#     Generate nodes and edges from the DataFrame. The user can select a column to use as the edge label.
#     Color specific nodes based on CLUSTERID1 and CLUSTERID2 from filtered_df.
#     Parameters:
#     df (pd.DataFrame): The DataFrame containing the node and edge data.
#     filtered_df (pd.DataFrame): The DataFrame containing a single row with CLUSTERID1 and CLUSTERID2.
#     edge_label_column (str): The name of the column from df to use as the label for the edges.
    
#     Returns:
#     nodes, edges: Lists of nodes and edges for the graph.
#     """
#     nodes = []
#     edges = []
#     added_nodes = set()  # Set to track added nodes to avoid duplicates

#     # Get CLUSTERID1 and CLUSTERID2 from filtered_df
#     id1 = str(filtered_df['CLUSTERID1'].iloc[0])
#     id2 = str(filtered_df['CLUSTERID2'].iloc[0])

#     # Add the blue node for CLUSTERID1 if it hasn't been added
#     if id1 not in added_nodes:
#         nodes.append(Node(id=id1,
#                           label="Source",  # Adjust label if needed
#                           size=20, 
#                           color="blue",
#                           ))
#         added_nodes.add(id1)

#     # Add the red node for CLUSTERID2 if it hasn't been added
#     if id2 not in added_nodes:
#         nodes.append(Node(id=id2,
#                           label="Target",  # Adjust label if needed
#                           size=20, 
#                           color="red",
#                           ))
#         added_nodes.add(id2)

#     for _, row in df.iterrows():
#         clusterid1 = str(row['CLUSTERID1'])
#         clusterid2 = str(row['CLUSTERID2'])
#         abs_chemprop = row['abs_ChemProp2']
#         sign_chemprop2 = row['Sign_ChemProp2']
#         deltamz = row['DeltaMZ']
#         id1_name = row['ID1_name'] if 'ID1_name' in row else clusterid1
#         id2_name = row['ID2_name'] if 'ID2_name' in row else clusterid2
#         mz1 = row['ID1_mz']
#         mz2 = row['ID2_mz']
#         rt1 = row['ID1_RT']
#         rt2 = row['ID2_RT']

#         # Get the edge label from the user-selected column
#         edge_label_value = row[edge_label_column]

#         # Color node1 and node2 based on id1 and id2
#         color_1 = "blue" if clusterid1 == id1 else "lightgray"
#         color_2 = "lightgray" if clusterid2 != id2 else "red"

#         # color_1 = "blue" if clusterid1 == id1 else "lightgray"
#         # color_2 = "red" if clusterid2 == id2 else "lightgray"

#         # Add nodes if not already added
#         if clusterid1 not in added_nodes:
#             nodes.append(Node(id=clusterid1, 
#                               label=f"{mz1:.2f}", 
#                               size=20, 
#                               color=color_1,
#                               title=f"ID: {clusterid1}\n Name: {id1_name}\n m/z: {mz1}\n RT: {rt1}"))
#             added_nodes.add(clusterid1)
        
#         if clusterid2 not in added_nodes:
#             nodes.append(Node(id=clusterid2,
#                               label=f"{mz2:.2f}",
#                               size=20, 
#                               color=color_2,
#                               title=f"ID: {clusterid2}\nName: {id2_name}\nm/z: {mz2}\nRT: {rt2}"))
#             added_nodes.add(clusterid2)

#         # Add edge with arrow based on abs_chemprop and sign_chemprop2
#         if abs_chemprop > 0:
#            if sign_chemprop2 == 1:
#                edges.append(Edge(source=clusterid1, 
#                                  target=clusterid2,label=f"{edge_label_value:.2f}", 
#                                  color="orange", 
#                                  arrow=True, 
#                                  font={"size": 10})
#                                  )
                
#            elif sign_chemprop2 == -1:
#                edges.append(Edge(source=clusterid2, 
#                                  target=clusterid1, 
#                                  label=f"{edge_label_value:.2f}", 
#                                  color="orange", arrow=True, 
#                                  font={"size": 10})
#                                  )
    
#     # Update the red node with the correct label and title information
#     for node in nodes:
#         if node.id == id1:
#             # Find the row in df with the correct information for id1
#             row = df[(df['CLUSTERID1'] == int(id1))].iloc[0]
#             node.label = f"{row['ID1_mz']:.2f}"
#             node.title = f"ID: {id1}\n Name: {row['ID1_name']}\n m/z: {row['ID1_mz']}\n RT: {row['ID1_RT']}"
#         if node.id == id2:
#             # Find the row in df with the correct information for id2
#             row = df[(df['CLUSTERID2'] == int(id2))].iloc[0]
#             node.label = f"{row['ID2_mz']:.2f}"
#             node.title = f"ID: {id2}\n Name: {row['ID2_name']}\n m/z: {row['ID2_mz']}\n RT: {row['ID2_RT']}"
#             break

#     return nodes, edges