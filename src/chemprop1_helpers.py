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
import plotly.express as px
import networkx as nx
from streamlit_agraph import agraph, Node, Edge, Config



def subset_levels(metatable, feature_table):
    if metatable.empty or feature_table.empty:
        st.warning("Metadata or Feature table is empty. Please check the data before subsetting.")
        return metatable, feature_table

    # User selects the attributes to subset
    all_columns = metatable.columns.tolist()
    selected_attributes = st.multiselect("Select multiple attributes to subset:", all_columns, default=all_columns[0])

    for attribute in selected_attributes:
        levels_cdtn = metatable[attribute].unique()  # unique levels for the selected attribute

        # User decision to keep or exclude conditions
        c1,c2 = st.columns(2)
        read_cdtn = c1.radio(f"Do you want to keep or exclude conditions for '{attribute}'?", 
                             ('Keep', 'Exclude'), 
                             key = f'radio_{attribute}', )

        # User selects conditions to keep or exclude
        selected_conditions = c2.multiselect(f"Select conditions to {read_cdtn.lower()} for '{attribute}':", 
                                             levels_cdtn, 
                                             key = f'multiselect_{attribute}')

        if len(selected_conditions) > 0:
             # Filter the metatable based on the user's selection
            if read_cdtn == 'Keep':
                metatable = metatable[metatable[attribute].isin(selected_conditions)]
            elif read_cdtn == 'Exclude':
                metatable = metatable[~metatable[attribute].isin(selected_conditions)]

    # Check if subset_metadata index is in feature table columns
    if set(metatable.index).issubset(set(feature_table.columns)):
        # Subset the feature table based on the subset metadata index
        subset_feature_table = feature_table[metatable.index]

        with st.expander(f"Metadata after subsetting {metatable.shape}"):
            st.dataframe(metatable)
        
        with st.expander(f"Feature table after subsetting {subset_feature_table.shape}"):
            st.dataframe(subset_feature_table)
    
    else:
        # Handle the case where some metadata indices are not in feature table columns
        missing_columns = set(metatable.index) - set(feature_table.columns)
        st.warning(f"Some metadata indices are not present as columns in the feature table: {missing_columns}")
        subset_feature_table = feature_table.copy()

    # Return the subsetted metadata and feature table
    return metatable, subset_feature_table


# Function to check if an attribute contains numerical values
def contains_convertible_values(series):
    for value in series.unique():
        # Try to convert each unique value to a number
        numeric_part = re.sub(r"[^\d.]", "", str(value))
        numeric_value = float(numeric_part) if numeric_part else 0.1

         # If the numeric value is zero, change it to 0.1
        if numeric_value == 0:
            numeric_value = 0.1
        try:
            float(numeric_part)
            return True  # Return True if conversion is successful
        except ValueError:
            continue  # Skip if conversion fails
    return False  # Return False if no values could be converted


def convert_series_values(series):
    converted_series = series.copy()  # Copy the series to avoid modifying the original one
    for index, value in series.items():
        # Strip away non-numeric characters
        numeric_part = re.sub(r"[^\d.]", "", str(value))

        # Convert to float; if conversion fails or results in zero, set to 0.0
        try:
            numeric_value = float(numeric_part)
            if numeric_value == 0:
                numeric_value = 0.0
        except ValueError:
            numeric_value = 0.0

        # Update the value in the copied series
        converted_series.at[index] = numeric_value

    return converted_series


# ChemProp analysis function
def run_chemprop1_analysis(desired_network_table, desired_feature_table, chemprop_metadata):
    # Initialize ChemProp1 score column with NA values (np.nan)
    ChemProp1 = np.full(len(desired_network_table), np.nan)

    for i in range(len(desired_network_table)):
        try:
            names_nodes = list(desired_feature_table.index)
            feature_ids = [names_nodes.index(desired_network_table.iloc[i, 0]), 
                           names_nodes.index(desired_network_table.iloc[i, 1])]
            
            # Subset feature data
            x = desired_feature_table.iloc[feature_ids, :]

            # Get metadata timepoints
            chemprop_sample_names = chemprop_metadata.index
            chemprop_timepoint = chemprop_metadata.iloc[:, 0]
            
            # Reorder `x` based on sample names in metadata
            reorder_id = [list(desired_feature_table.columns).index(name) for name in chemprop_sample_names]
            reordered_x = x.iloc[:, reorder_id].T
            reordered_x = pd.concat([
                chemprop_timepoint.reset_index(drop=True), 
                reordered_x.reset_index(drop=True)
                ], axis=1)
            
            # Sort `reordered_x` by the first column (0th index) in ascending order
            reordered_x = reordered_x.sort_values(by=reordered_x.columns[0], ascending=True).reset_index(drop=True)

            # Subset features by cluster IDs
            Feature1 = reordered_x.iloc[:,1]
            Feature2 = reordered_x.iloc[:,2]
            timepoints = reordered_x.iloc[:,0].unique()
            timepoint1, timepoint2 = timepoints[0], timepoints[1]

            reordered_x = reordered_x.rename(columns={
                reordered_x.columns[0]: 'Timepoint',  # Rename the first column to 'Timepoint'
                reordered_x.columns[1]: 'Feature1',   # Rename the second column to 'Feature1'
                reordered_x.columns[2]: 'Feature2'    # Rename the third column to 'Feature2'
            })

            # Subset values based on timepoints
            A1 = reordered_x.loc[reordered_x['Timepoint'] == timepoint1, 'Feature1']
            A2 = reordered_x.loc[reordered_x['Timepoint'] == timepoint2, 'Feature1']
            B1 = reordered_x.loc[reordered_x['Timepoint'] == timepoint1, 'Feature2']
            B2 = reordered_x.loc[reordered_x['Timepoint'] == timepoint2, 'Feature2']

            # Provide a constant k to avoid division by zero
            if (
                st.session_state['blank_removal_done'] == True 
                and st.session_state['imputation_done'] == False 
                and st.session_state['normalization_done'] == False
            ):
                k = 1e-10 #1.0 e−10
            elif (
                st.session_state['imputation_done'] == True 
                or st.session_state['normalization_done'] == True
            ):
                k = 0

            # Calculate score with log ratio
            top = ((A1 + k) / (B1 + k)).mean()
            bottom = (((A2 + k) / (B2 + k))).mean()
            ChemProp1[i] = np.log(top/bottom)

        except ValueError:
            # If the feature ID is not found, set ChemProp1 to NaN
            ChemProp1[i] = np.nan

    # Combine results into a new DataFrame
    new_network_table = desired_network_table.copy()
    new_network_table['ChemProp1'] = ChemProp1

    # Add absolute value columns
    abs_values = new_network_table[['ChemProp1']].abs()
    abs_values.columns = [f"abs_{col}" for col in abs_values.columns]
    abs_values = abs_values.apply(pd.to_numeric, errors='coerce')

    # Add sign of ChemProp1
    new_network_table['Sign_ChemProp1'] = np.sign(new_network_table['ChemProp1'])
    
    # Concatenate absolute values
    ChemProp1_file = pd.concat([new_network_table, abs_values], axis=1)

    return ChemProp1_file


##################

# Function to strip non-numeric characters
def strip_non_numeric(val):
    if isinstance(val, str):  # Only apply the regex to string values
        return re.sub(r'[^\d.]+', '', val)  # Remove everything except digits and dots
    return val  # If it's not a string (e.g., already numeric), return it as-is


def add_names_to_chemprop(edge_df, gnps_df):
    """
    Adds 'ID1_name', 'ID2_name', 'ID1_mz', 'ID2_mz', 'ID1_RT', and 'ID2_RT' columns to the edge_df
    by matching 'CLUSTERID1' and 'CLUSTERID2' with 'cluster index' in gnps_df and adding the 
    corresponding 'Compound Name', 'parent mass', and 'RTMean' to the new columns.

    Args:
        edge_df (pd.DataFrame): The ChemProp2_scores dataframe.
        gnps_df (pd.DataFrame): The an_gnps dataframe.

    Returns:
        pd.DataFrame: The updated dataframe with the new columns added.
    """
    # Determine the column renaming based on the existing columns in gnps_df
    if {'cluster index', 'parent mass', 'RTMean'}.issubset(gnps_df.columns):
        gnps_df = gnps_df.rename(columns={'cluster index': 'CLUSTERID', 'parent mass': 'mz', 'RTMean': 'RT'})
    elif {'#Scan#', 'Precursor_MZ', 'RT_Query'}.issubset(gnps_df.columns):
        gnps_df = gnps_df.rename(columns={'#Scan#': 'CLUSTERID', 'Precursor_MZ': 'mz', 'RT_Query': 'RT'})
        # 'RT' is not present, so we drop it if not available
        if 'RT' in gnps_df.columns:
            gnps_df = gnps_df[['CLUSTERID', 'Compound_Name', 'mz', 'RT']]
        else:
            gnps_df = gnps_df[['CLUSTERID', 'Compound_Name', 'mz']]

    # Verify required columns are in gnps_df after renaming
    if 'Compound_Name' not in gnps_df.columns:
        raise KeyError("The 'Compound_Name' column is missing from the annotation file. Check the input data.")
    
    # Initialize the columns for mz and RT in edge_df based on gnps_df availability
    cols_to_merge = ['CLUSTERID', 'Compound_Name']
    if 'mz' in gnps_df.columns:
        cols_to_merge.append('mz')
    if 'RT' in gnps_df.columns:
        cols_to_merge.append('RT')

    # Merge for CLUSTERID1
    merged_df = edge_df.merge(gnps_df[cols_to_merge],
                              left_on='CLUSTERID1',
                              right_on='CLUSTERID',
                              how='left')
    # Rename columns from gnps_df to indicate they correspond to CLUSTERID1
    rename_columns = {'Compound_Name': 'ID1_name'}
    if 'mz' in gnps_df.columns:
        rename_columns['mz'] = 'ID1_mz'
    if 'RT' in gnps_df.columns:
        rename_columns['RT'] = 'ID1_RT'
    merged_df = merged_df.rename(columns=rename_columns)
    merged_df = merged_df.drop(columns=['CLUSTERID'])

    # Merge for CLUSTERID2
    merged_df = merged_df.merge(gnps_df[cols_to_merge],
                                left_on='CLUSTERID2',
                                right_on='CLUSTERID',
                                how='left')
    # Rename columns from gnps_df to indicate they correspond to CLUSTERID2
    rename_columns = {'Compound_Name': 'ID2_name'}
    if 'mz' in gnps_df.columns:
        rename_columns['mz'] = 'ID2_mz'
    if 'RT' in gnps_df.columns:
        rename_columns['RT'] = 'ID2_RT'
    merged_df = merged_df.rename(columns=rename_columns)

    # Drop the extra 'CLUSTERID' column generated during merging
    merged_df = merged_df.drop(columns=['CLUSTERID'], errors='ignore')

    # If RT was not present in gnps_df, add 'ID1_RT' and 'ID2_RT' columns with zeros
    if 'ID1_RT' not in merged_df.columns:
        merged_df['ID1_RT'] = 0
    if 'ID2_RT' not in merged_df.columns:
        merged_df['ID2_RT'] = 0

    return merged_df



def generate_graphml_with_secondary_edges_chemprop1(df, output_file):
    """
    Generate a GraphML file from the ChemProp2_scores DataFrame
    
    Parameters:
    df (pd.DataFrame): The ChemProp2_scores DataFrame containing node and edge data.
    output_file (str): The name of the GraphML file to generate.
    """
    # Create a directed graph allowing multiple edges
    G = nx.MultiDiGraph()

    # Check for 'chemprop_ft' in session state for additional node data
    node_table = st.session_state['chemprop1_ft'] if 'chemprop1_ft' in st.session_state else None

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
            'ChemProp1': row['ChemProp1'],
            'label': f"∆mz: {row['DeltaMZ']:.2f}",
            'color': "black"
        }

        # Add primary edge
        G.add_edge(clusterid1, clusterid2, key="primary", **primary_edge_attributes)

        # Define secondary edge attributes based on abs_ChemProp2 and Sign_ChemProp2
        abs_chemprop1 = row['abs_ChemProp1']
        sign_chemprop1 = row['Sign_ChemProp1']
        secondary_edge_attributes = {
                'weight': abs_chemprop1,
                'label': f"{abs_chemprop1:.2f}",
                'color': "red",
                'abs_ChemProp1': row['abs_ChemProp1']
                }
        # Use Sign_ChemProp2 to determine direction of the secondary edge
        if sign_chemprop1 == 1:
            G.add_edge(clusterid1, clusterid2, key="secondary", **secondary_edge_attributes, arrow=True)
        elif sign_chemprop1 == -1:
            G.add_edge(clusterid2, clusterid1, key="secondary", **secondary_edge_attributes, arrow=True)

    # Write the graph to GraphML
    nx.write_graphml(G, output_file)

    return output_file


def generate_graphml_zip_chemprop1():
    if 'ChemProp1_scores' in st.session_state:
        df = st.session_state['ChemProp1_scores'].copy()

        chemprop1_graphml_filename = "chemprop1_graph.graphml"
        style_filename = "resources/ChemProp1_styles.xml"

        with tempfile.TemporaryDirectory() as temp_dir:
            chemprop1_graphml_path = os.path.join(temp_dir, chemprop1_graphml_filename)
            generate_graphml_with_secondary_edges_chemprop1(df, chemprop1_graphml_path)

            # Full path to styles.xml
            style_path = os.path.join(os.path.dirname(__file__), style_filename)
            
            # Check if the styles.xml file exists at the calculated path
            if not os.path.exists(style_path):
                st.error("The styles.xml file for ChemProp1 is not found in the 'resources' folder.")
                return

            # Create the zip file in the temporary directory
            zip_path = os.path.join(temp_dir, 'chemprop1_graph_with_styles.zip')
            with zipfile.ZipFile(zip_path, 'w') as zipf:
                zipf.write(chemprop1_graphml_path, arcname=chemprop1_graphml_filename)
                zipf.write(style_path, arcname="styles.xml")

            # Provide the zip file for download without leaving it in the repo
            with open(zip_path, 'rb') as f:
                st.download_button(
                    label="Download GraphML and Styles as ZIP",
                    data=f,
                    file_name='chemprop1_graph_with_styles.zip',
                    mime='application/zip'
                )


def plot_intensity_trends_single_row_chemprop1(row, feature_table, metadata):
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
    chem_score = row['ChemProp1']

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
        long_df['Timepoint'] = long_df['Timepoint'].astype(str)

        # Create box plot without mean lines
        fig = px.box(
            long_df,
            x='Timepoint',  # X-axis for time points
            y='Value',  # Y-axis for the intensity values of both features
            color='Feature',  # Color by 'Feature' to distinguish Feature1 and Feature2
            template="plotly_white",
            title=f"Intensity Trend Box Plot {clusterid1} vs {clusterid2}",
            labels={'Timepoint': 'Time', 'Value': 'Feature Intensity'},
            width=800,
            height=600,
            points="all",
        )

        # Add subtitle for ChemProp2 score using annotation
        fig.add_annotation(
           text=f"ChemProp1 score = {chem_score:.2f}<br>{name_id1} <b>VS</b> {name_id2}",
           xref="paper", yref="paper",
           x=0.5, y=1.1,  # Positioning the subtitle above the plot
           showarrow=False,
           font=dict(size=12)
        )

        # Return the Plotly figure object
        return fig

def generate_graph_from_df_chemprop1(df, filtered_df, edge_label_column):
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

    # Add the blue node for CLUSTERID1 if it hasn't been added
    if id1 not in added_nodes:
        nodes.append(Node(id=id1,
                          label="Source",  # Adjust label if needed
                          size=20, 
                          color="blue",
                          ))
        added_nodes.add(id1)

    # Add the red node for CLUSTERID2 if it hasn't been added
    if id2 not in added_nodes:
        nodes.append(Node(id=id2,
                          label="Target",  # Adjust label if needed
                          size=20, 
                          color="red",
                          ))
        added_nodes.add(id2)

    for _, row in df.iterrows():
        clusterid1 = str(row['CLUSTERID1'])
        clusterid2 = str(row['CLUSTERID2'])
        abs_chemprop = row['abs_ChemProp1']
        sign_chemprop = row['Sign_ChemProp1']
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
        color_2 = "lightgray" if clusterid2 != id2 else "red"

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
                              title=f"ID: {clusterid2}\n Name: {id2_name}\n m/z: {mz2}\n RT: {rt2}"))
            added_nodes.add(clusterid2)

        # Add edge with arrow based on abs_chemprop and sign_chemprop2
        if abs_chemprop > 0:
           if sign_chemprop == 1:
               edges.append(Edge(source=clusterid1, 
                                 target=clusterid2,
                                 label=f"{edge_label_value:.2f}", 
                                 color="orange", 
                                 arrow=True, 
                                 font={"size": 10})
                                 )
                
           elif sign_chemprop == -1:
               edges.append(Edge(source=clusterid2, 
                                 target=clusterid1, 
                                 label=f"{edge_label_value:.2f}", 
                                 color="orange", 
                                 arrow=True, 
                                 font={"size": 10})
                                 )

    # Update the red node with the correct label and title information
    for node in nodes:
        if node.id == id1:
            # Find the row in df with the correct information for id1
            row = df[(df['CLUSTERID1'] == int(id1))].iloc[0]
            node.label = f"{row['ID1_mz']:.2f}"
            node.title = f"ID: {id1}\n Name: {row['ID1_name']}\n m/z: {row['ID1_mz']}\n RT: {row['ID1_RT']}"
        if node.id == id2:
            # Find the row in df with the correct information for id2
            row = df[(df['CLUSTERID2'] == int(id2))].iloc[0]
            node.label = f"{row['ID2_mz']:.2f}"
            node.title = f"ID: {id2}\n Name: {row['ID2_name']}\n m/z: {row['ID2_mz']}\n RT: {row['ID2_RT']}"
            break

    return nodes, edges


# def plot_intensity_trends_single_row_chemprop1(row, feature_table, metadata):
#     """
#     Plots intensity trends for a single row containing CLUSTERID1 and CLUSTERID2,
#     where the x-axis represents unique timepoints from the metadata.

#     Parameters:
#     row (pd.Series): A single row from the DataFrame containing CLUSTERID1 and CLUSTERID2.
#     feature_table (pd.DataFrame): The feature table containing intensity values.
#     metadata (pd.DataFrame): The metadata containing timepoints (single column).
#     """
#     # Extract cluster IDs and relevant information from the row
#     clusterid1 = row['CLUSTERID1']
#     clusterid2 = row['CLUSTERID2']
#     name_id1 = row['ID1_name']
#     name_id2 = row['ID2_name']
#     chem_score = row['ChemProp1']

#     # Extract corresponding rows from feature_table for CLUSTERID1 and CLUSTERID2
#     if clusterid1 in feature_table.index and clusterid2 in feature_table.index:
#         intensity_1 = feature_table.loc[clusterid1]
#         intensity_2 = feature_table.loc[clusterid2]

#         # Align intensity data with the timepoints from the metadata
#         intensity_1_with_timepoints = pd.concat([intensity_1, metadata], axis=1)
#         intensity_2_with_timepoints = pd.concat([intensity_2, metadata], axis=1)

#         # Merge intensities based on timepoints
#         merged_total = pd.merge(intensity_1_with_timepoints,
#                                 intensity_2_with_timepoints,
#                                 left_index=True,
#                                 right_index=True,
#                                 suffixes=('_CLUSTERID1', '_CLUSTERID2'))

#         # Ensure all column names are strings to safely apply 'endswith'
#         merged_total.columns = merged_total.columns.map(str)

#         # Dropping columns that end with '_CLUSTERID1' from the dataframe
#         merged_total = merged_total.drop(columns=[col for col in merged_total.columns if col.endswith('_CLUSTERID1')])

#         # Group by timepoints and calculate the average intensity for each timepoint
#         intensity_1_grouped = intensity_1_with_timepoints.groupby(metadata.columns[0]).mean()
#         intensity_2_grouped = intensity_2_with_timepoints.groupby(metadata.columns[0]).mean()

#         # Merge the two grouped DataFrames by index (timepoints)
#         mean_intensities = pd.merge(intensity_1_grouped, 
#                                     intensity_2_grouped, 
#                                     left_index=True, 
#                                     right_index=True, 
#                                     suffixes=('_CLUSTERID1', '_CLUSTERID2'))

#         # Reset the index if 'Timepoint' is currently the row index
#         mean_intensities = mean_intensities.reset_index()

#         # Ensure there are exactly 3 columns to rename
#         if mean_intensities.shape[1] == 3:
#             mean_intensities = mean_intensities.rename(columns={
#                 mean_intensities.columns[0]: 'Timepoint',  # Rename the first column to 'Timepoint'
#                 mean_intensities.columns[1]: 'Feature1',   # Rename the second column to 'Feature1'
#                 mean_intensities.columns[2]: 'Feature2'    # Rename the third column to 'Feature2'
#             })

#         # Convert data to long format for plotting
#         long_df = pd.melt(merged_total, 
#                           id_vars=[merged_total.columns[2]], 
#                           value_vars=[merged_total.columns[0], 
#                                       merged_total.columns[1]], 
#                                       var_name='Feature', 
#                                       value_name='Value'
#                           )
        
#         long_df['Feature'] = 'Feature_' + long_df['Feature'].astype(str)
#         long_df.rename(columns={long_df.columns[0]: 'Timepoint'}, inplace=True)

#         # Create scatter plot with lines connecting the points
#         fig = px.scatter(
#             long_df,
#             x='Timepoint',  # X-axis for time points
#             y='Value',  # Y-axis for the intensity values of both features
#             color='Feature',  # Color by 'Feature' to distinguish Feature1 and Feature2
#             template="plotly_white",
#             title=f"Intensity Trend Plot {clusterid1} vs {clusterid2}",
#             labels={'Timepoint': 'Time', 'Value': 'Feature Intensity'},  # Ensure 'Timepoint' label matches your DataFrame
#             width=800,
#             height=600,
#         )

#         # Add mean lines for Feature1 and Feature2
#         fig.add_scatter(
#             x=mean_intensities['Timepoint'], 
#             y=mean_intensities['Feature1'], 
#             mode='lines', 
#             name=f"Mean {clusterid1}",  # Label for the mean line of Feature1
#             line=dict(color='blue', dash='dash')  # Customize color and style of the line
#         )

#         fig.add_scatter(
#             x=mean_intensities['Timepoint'], 
#             y=mean_intensities['Feature2'], 
#             mode='lines', 
#             name=f"Mean {clusterid2}",  # Label for the mean line of Feature2
#             line=dict(color='red', dash='dash')  # Customize color and style of the line
#         )

#         # Add subtitle for ChemProp2 score using annotation
#         fig.add_annotation(
#            text=f"ChemProp1 score = {chem_score:.2f}<br>{name_id1} <b>VS</b> {name_id2}",
#            xref="paper", yref="paper",
#            x=0.5, y=1.1,  # Positioning the subtitle above the plot
#            showarrow=False,
#            font=dict(size=12)
#         )

#         # Return the Plotly figure object
#         return fig

