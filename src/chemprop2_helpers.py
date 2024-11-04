import streamlit as st
import os
import zipfile
import tempfile
from io import BytesIO
from .common import *  # Importing common functionalities from the 'common' module
from src.cleanup import *
import pandas as pd
import numpy as np
import re
from scipy.stats import spearmanr, pearsonr
import plotly.graph_objects as go
import networkx as nx
from streamlit_agraph import agraph, Node, Edge, Config

def ChemProp2_scoring(desired_network_table, desired_feature_table, names_nodes, chemprop_metadata):
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
    
    # Add sign of ChemProp2
    new_network_table['Sign_ChemProp2'] = np.sign(new_network_table['ChemProp2'])
    
    # Concatenate absolute values
    ChemProp2_file = pd.concat([new_network_table, abs_values], axis=1)

    return ChemProp2_file

def plot_intensity_trends_single_row_chemprop2(row, feature_table, metadata):
    """
    Plots intensity trends for a single row containing CLUSTERID1 and CLUSTERID2,
    where the x-axis represents unique timepoints from the metadata.

    Parameters:
    row (pd.Series): A single row from the DataFrame containing CLUSTERID1 and CLUSTERID2.
    feature_table (pd.DataFrame): The feature table containing intensity values.
    metadata (pd.DataFrame): The metadata containing timepoints (single column).
    """
    # Extract cluster IDs and relevant information from the row
    clusterid1 = row['CLUSTERID1']
    clusterid2 = row['CLUSTERID2']
    name_id1 = row['ID1_name']
    name_id2 = row['ID2_name']
    chem_score = row['ChemProp2']

    # Extract corresponding rows from feature_table for CLUSTERID1 and CLUSTERID2
    if clusterid1 in feature_table.index and clusterid2 in feature_table.index:
        intensity_1 = feature_table.loc[clusterid1]
        intensity_2 = feature_table.loc[clusterid2]

        # Align intensity data with the timepoints from the metadata
        intensity_1_with_timepoints = pd.concat([intensity_1, metadata], axis=1)
        intensity_2_with_timepoints = pd.concat([intensity_2, metadata], axis=1)

        # Merge intensities based on timepoints
        merged_total = pd.merge(intensity_1_with_timepoints,
                                intensity_2_with_timepoints,
                                left_index=True,
                                right_index=True,
                                suffixes=('_CLUSTERID1', '_CLUSTERID2'))

        # Ensure all column names are strings to safely apply 'endswith'
        merged_total.columns = merged_total.columns.map(str)

        # Dropping columns that end with '_CLUSTERID1' from the dataframe
        merged_total = merged_total.drop(columns=[col for col in merged_total.columns if col.endswith('_CLUSTERID1')])

        # Group by timepoints and calculate the average intensity for each timepoint
        intensity_1_grouped = intensity_1_with_timepoints.groupby(metadata.columns[0]).mean()
        intensity_2_grouped = intensity_2_with_timepoints.groupby(metadata.columns[0]).mean()

        # Merge the two grouped DataFrames by index (timepoints)
        mean_intensities = pd.merge(intensity_1_grouped, 
                                    intensity_2_grouped, 
                                    left_index=True, 
                                    right_index=True, 
                                    suffixes=('_CLUSTERID1', '_CLUSTERID2'))

        # Reset the index if 'Timepoint' is currently the row index
        mean_intensities = mean_intensities.reset_index()

        # Ensure there are exactly 3 columns to rename
        if mean_intensities.shape[1] == 3:
            mean_intensities = mean_intensities.rename(columns={
                mean_intensities.columns[0]: 'Timepoint',  # Rename the first column to 'Timepoint'
                mean_intensities.columns[1]: 'Feature1',   # Rename the second column to 'Feature1'
                mean_intensities.columns[2]: 'Feature2'    # Rename the third column to 'Feature2'
            })

        # Convert data to long format for plotting
        long_df = pd.melt(merged_total, 
                          id_vars=[merged_total.columns[2]], 
                          value_vars=[merged_total.columns[0], 
                                      merged_total.columns[1]], 
                                      var_name='Feature', 
                                      value_name='Value'
                          )
        
        long_df['Feature'] = 'Feature_' + long_df['Feature'].astype(str)
        long_df.rename(columns={long_df.columns[0]: 'Timepoint'}, inplace=True)

        # Create scatter plot with lines connecting the points
        fig = px.scatter(
            long_df,
            x='Timepoint',  # X-axis for time points
            y='Value',  # Y-axis for the intensity values of both features
            color='Feature',  # Color by 'Feature' to distinguish Feature1 and Feature2
            template="plotly_white",
            title=f"Intensity Trend Plot {clusterid1} vs {clusterid2}",
            labels={'Timepoint': 'Time', 'Value': 'Feature Intensity'},  # Ensure 'Timepoint' label matches your DataFrame
            width=800,
            height=600,
        )

        # Add mean lines for Feature1 and Feature2
        fig.add_scatter(
            x=mean_intensities['Timepoint'], 
            y=mean_intensities['Feature1'], 
            mode='lines', 
            name=f"Mean {clusterid1}",  # Label for the mean line of Feature1
            line=dict(color='blue', dash='dash')  # Customize color and style of the line
        )

        fig.add_scatter(
            x=mean_intensities['Timepoint'], 
            y=mean_intensities['Feature2'], 
            mode='lines', 
            name=f"Mean {clusterid2}",  # Label for the mean line of Feature2
            line=dict(color='red', dash='dash')  # Customize color and style of the line
        )

        # Add subtitle for ChemProp2 score using annotation
        fig.add_annotation(
           text=f"ChemProp2 score = {chem_score:.2f}<br>{name_id1} <b>VS</b> {name_id2}",
           xref="paper", yref="paper",
           x=0.5, y=1.1,  # Positioning the subtitle above the plot
           showarrow=False,
           font=dict(size=12)
        )

        # Return the Plotly figure object
        return fig


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


def generate_graph_from_df_chemprop2(df, filtered_df, edge_label_column):
    """
    Generate nodes and edges from the DataFrame. The user can select a column to use as the edge label.
    Color specific nodes based on CLUSTERID1 and CLUSTERID2 from filtered_df.

    Parameters:
    df (pd.DataFrame): The DataFrame containing the node and edge data.
    filtered_df (pd.DataFrame): The DataFrame containing a single row with CLUSTERID1 and CLUSTERID2.
    edge_label_column (str): The name of the column from df to use as the label for the edges.
    
    Returns:
    nodes, edges: Lists of nodes and edges for the graph.
    """
    nodes = []
    edges = []
    added_nodes = set()  # Set to track added nodes to avoid duplicates

    # Get CLUSTERID1 and CLUSTERID2 from filtered_df
    id1 = str(filtered_df['CLUSTERID1'].iloc[0])
    id2 = str(filtered_df['CLUSTERID2'].iloc[0])

    for _, row in df.iterrows():
        clusterid1 = str(row['CLUSTERID1'])
        clusterid2 = str(row['CLUSTERID2'])
        abs_chemprop = row['abs_ChemProp2']
        sign_chemprop2 = row['Sign_ChemProp2']
        deltamz = row['DeltaMZ']
        id1_name = row['ID1_name'] if 'ID1_name' in row else clusterid1
        id2_name = row['ID2_name'] if 'ID2_name' in row else clusterid2
        mz1 = row['ID1_mz']
        mz2 = row['ID2_mz']
        rt1 = row['ID1_RT']
        rt2 = row['ID2_RT']

        # Get the edge label from the user-selected column
        edge_label_value = row[edge_label_column]

        # Color node1 and node2 based on id1 and id2
        color_1 = "blue" if clusterid1 == id1 else "lightgray"
        color_2 = "red" if clusterid2 == id2 else "lightgray"

        # Add nodes if not already added
        if clusterid1 not in added_nodes:
            nodes.append(Node(id=clusterid1, 
                              label=f"{mz1:.2f}", 
                              size=20, 
                              color=color_1,
                              title=f"ID: {clusterid1}\n Name: {id1_name}\n m/z: {mz1}\n RT: {rt1}"))
            added_nodes.add(clusterid1)
        
        if clusterid2 not in added_nodes:
            nodes.append(Node(id=clusterid2,
                              label=f"{mz2:.2f}",
                              size=20, 
                              color=color_2,
                              title=f"ID: {clusterid2}\nName: {id2_name}\nm/z: {mz2}\nRT: {rt2}"))
            added_nodes.add(clusterid2)
        
        # Add edge with arrow based on abs_chemprop and sign_chemprop2
        if abs_chemprop > 0:
           if sign_chemprop2 == 1:
               edges.append(Edge(source=clusterid1, target=clusterid2,label=f"{edge_label_value:.2f}", color="orange", arrow=True, font={"size": 10}))
                
           elif sign_chemprop2 == -1:
               edges.append(Edge(source=clusterid2, target=clusterid1, label=f"{edge_label_value:.2f}", color="orange", arrow=True, font={"size": 10}))

    return nodes, edges
