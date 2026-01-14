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
import urllib.parse as up


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


# def add_names_to_chemprop(edge_df, gnps_df):
#     """
#     Adds 'ID1_name', 'ID2_name', 'ID1_mz', 'ID2_mz', 'ID1_RT', and 'ID2_RT' columns to the edge_df
#     by matching 'CLUSTERID1' and 'CLUSTERID2' with 'cluster index' in gnps_df and adding the 
#     corresponding 'Compound Name', 'parent mass', and 'RTMean' to the new columns.

#     Args:
#         edge_df (pd.DataFrame): The ChemProp2_scores dataframe.
#         gnps_df (pd.DataFrame): The an_gnps dataframe.

#     Returns:
#         pd.DataFrame: The updated dataframe with the new columns added.
#     """
#     # Determine the column renaming based on the existing columns in gnps_df
#     if {'cluster index', 'parent mass', 'RTMean'}.issubset(gnps_df.columns):
#         gnps_df = gnps_df.rename(columns={'cluster index': 'CLUSTERID', 'parent mass': 'mz', 'RTMean': 'RT'})
#     elif {'#Scan#', 'Precursor_MZ', 'RT_Query'}.issubset(gnps_df.columns):
#         gnps_df = gnps_df.rename(columns={'#Scan#': 'CLUSTERID', 'Precursor_MZ': 'mz', 'RT_Query': 'RT'})
#         # 'RT' is not present, so we drop it if not available
#         if 'RT' in gnps_df.columns:
#             gnps_df = gnps_df[['CLUSTERID', 'Compound_Name', 'mz', 'RT']]
#         else:
#             gnps_df = gnps_df[['CLUSTERID', 'Compound_Name', 'mz']]

#     # Verify required columns are in gnps_df after renaming
#     if 'Compound_Name' not in gnps_df.columns:
#         raise KeyError("The 'Compound_Name' column is missing from the annotation file. Check the input data.")
    
#     # Initialize the columns for mz and RT in edge_df based on gnps_df availability
#     cols_to_merge = ['CLUSTERID', 'Compound_Name']
#     if 'mz' in gnps_df.columns:
#         cols_to_merge.append('mz')
#     if 'RT' in gnps_df.columns:
#         cols_to_merge.append('RT')

#     # Merge for CLUSTERID1
#     merged_df = edge_df.merge(gnps_df[cols_to_merge],
#                               left_on='CLUSTERID1',
#                               right_on='CLUSTERID',
#                               how='left')
#     # Rename columns from gnps_df to indicate they correspond to CLUSTERID1
#     rename_columns = {'Compound_Name': 'ID1_name'}
#     if 'mz' in gnps_df.columns:
#         rename_columns['mz'] = 'ID1_mz'
#     if 'RT' in gnps_df.columns:
#         rename_columns['RT'] = 'ID1_RT'
#     merged_df = merged_df.rename(columns=rename_columns)
#     merged_df = merged_df.drop(columns=['CLUSTERID'])

#     # Merge for CLUSTERID2
#     merged_df = merged_df.merge(gnps_df[cols_to_merge],
#                                 left_on='CLUSTERID2',
#                                 right_on='CLUSTERID',
#                                 how='left')
#     # Rename columns from gnps_df to indicate they correspond to CLUSTERID2
#     rename_columns = {'Compound_Name': 'ID2_name'}
#     if 'mz' in gnps_df.columns:
#         rename_columns['mz'] = 'ID2_mz'
#     if 'RT' in gnps_df.columns:
#         rename_columns['RT'] = 'ID2_RT'
#     merged_df = merged_df.rename(columns=rename_columns)

#     # Drop the extra 'CLUSTERID' column generated during merging
#     merged_df = merged_df.drop(columns=['CLUSTERID'], errors='ignore')

#     # If RT was not present in gnps_df, add 'ID1_RT' and 'ID2_RT' columns with zeros
#     if 'ID1_RT' not in merged_df.columns:
#         merged_df['ID1_RT'] = 0
#     if 'ID2_RT' not in merged_df.columns:
#         merged_df['ID2_RT'] = 0

#     return merged_df

def add_mz_rt_from_ft(
    edge_df: pd.DataFrame,
    ft_df: pd.DataFrame,
    *,
    id1_col: str = "CLUSTERID1",
    id2_col: str = "CLUSTERID2",
    ft_mz_col: str = "row m/z",
    ft_rt_col: str = "row retention time",
) -> pd.DataFrame:
    """
    Always add ID1_mz, ID1_RT, ID2_mz, ID2_RT from the feature table (ft_df).

    Assumes:
      - ft_df.index contains feature IDs (row ID)
      - ft_df has columns `row m/z` and `row retention time`
    """

    # --- sanity checks ---
    missing_edge = {id1_col, id2_col} - set(edge_df.columns)
    if missing_edge:
        raise KeyError(f"edge_df missing columns: {sorted(missing_edge)}")

    missing_ft = {ft_mz_col, ft_rt_col} - set(ft_df.columns)
    if missing_ft:
        raise KeyError(f"ft_df missing columns: {sorted(missing_ft)}")

    out = edge_df.copy()

    # --- normalize IDs (string join is safest) ---
    out[id1_col] = out[id1_col].astype("Int64").astype(str)
    out[id2_col] = out[id2_col].astype("Int64").astype(str)

    ft_meta = ft_df[[ft_mz_col, ft_rt_col]].copy()
    ft_meta.index = ft_meta.index.astype("Int64").astype(str)

    # --- merge ID1 ---
    out = out.merge(
        ft_meta.rename(columns={ft_mz_col: "ID1_mz", ft_rt_col: "ID1_RT"}),
        how="left",
        left_on=id1_col,
        right_index=True,
        validate="m:1",
    )

    # --- merge ID2 ---
    out = out.merge(
        ft_meta.rename(columns={ft_mz_col: "ID2_mz", ft_rt_col: "ID2_RT"}),
        how="left",
        left_on=id2_col,
        right_index=True,
        validate="m:1",
    )

    # --- add DeltaMZ if missing ---
    if "DeltaMZ" not in out.columns:
        out["DeltaMZ"] = (out["ID1_mz"] - out["ID2_mz"]).abs()

    return out

def add_names_from_gnps(
    edge_df: pd.DataFrame,
    gnps_df: pd.DataFrame,
    *,
    id1_col: str = "CLUSTERID1",
    id2_col: str = "CLUSTERID2",
) -> pd.DataFrame:
    """
    If gnps_df exists, add ID1_name and ID2_name from GNPS annotations.
    Only names are added here (no mz/RT).
    """
    if gnps_df is None or not isinstance(gnps_df, pd.DataFrame) or gnps_df.empty:
        return edge_df

    out = edge_df.copy()

    # Normalize GNPS schema -> columns: CLUSTERID, Compound_Name
    g = gnps_df.copy()
    if {"cluster index", "Compound_Name"}.issubset(g.columns):
        g = g.rename(columns={"cluster index": "CLUSTERID"})
    elif {"#Scan#", "Compound_Name"}.issubset(g.columns):
        g = g.rename(columns={"#Scan#": "CLUSTERID"})
    elif "CLUSTERID" not in g.columns:
        raise KeyError("gnps_df must contain 'cluster index' or '#Scan#' (or already 'CLUSTERID').")

    if "Compound_Name" not in g.columns:
        raise KeyError("gnps_df is missing 'Compound_Name' column.")

    g = g[["CLUSTERID", "Compound_Name"]].dropna(subset=["CLUSTERID"]).copy()

    # Normalize join keys
    out[id1_col] = out[id1_col].astype("Int64").astype(str)
    out[id2_col] = out[id2_col].astype("Int64").astype(str)
    g["CLUSTERID"] = g["CLUSTERID"].astype("Int64").astype(str)

    # Deduplicate GNPS by cluster id (keep first name)
    g = g.drop_duplicates(subset=["CLUSTERID"])

    # Merge ID1 name
    out = out.merge(
        g.rename(columns={"CLUSTERID": "_GNPS_ID", "Compound_Name": "ID1_name"}),
        how="left",
        left_on=id1_col,
        right_on="_GNPS_ID",
        validate="m:1",
    ).drop(columns=["_GNPS_ID"])

    # Merge ID2 name
    out = out.merge(
        g.rename(columns={"CLUSTERID": "_GNPS_ID", "Compound_Name": "ID2_name"}),
        how="left",
        left_on=id2_col,
        right_on="_GNPS_ID",
        validate="m:1",
    ).drop(columns=["_GNPS_ID"])

    return out


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

    # normalize index
    ft = feature_table.copy()
    ft.index = ft.index.astype(str)
    id1 = str(int(clusterid1))
    id2 = str(int(clusterid2))

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

    label1 = f"ID {clusterid1}"
    label2 = f"ID {clusterid2}"

    df = pd.DataFrame({
        "Timepoint": metadata.loc[common_samples, time_col].astype(str).values,
        label1: ft.loc[id1, common_samples].values,
        label2: ft.loc[id2, common_samples].values,
    })

    long_df = df.melt(
        id_vars="Timepoint",
        value_vars=[label1, label2],
        var_name="Feature",
        value_name="Value"
    )

    color_map = {
        label1: "blue",  # blue
        label2: "#d62728",  # red
    }

    fig = px.box(
        long_df,
        x="Timepoint",
        y="Value",
        color="Feature",
        color_discrete_map=color_map,
        template="plotly_white",
        title=f"Intensity Trend Box Plot {clusterid1} vs {clusterid2}",
        labels={"Timepoint": "Time", "Value": "Feature Intensity"},
        points="all",
    )

    score = row.get("ChemProp1", np.nan)
    name1 = row.get("ID1_name", "")
    name2 = row.get("ID2_name", "")

    fig.add_annotation(
        text=f"ChemProp1 score = {'NA' if pd.isna(score) else f'{score:.2f}'}<br>"
             f"{name1} <b>VS</b> {name2}",
        xref="paper", yref="paper",
        x=0.5, y=1.1,
        showarrow=False,
    )

    return fig, status



def generate_graph_from_df_chemprop1(df: pd.DataFrame, filtered_df: pd.DataFrame, edge_label_column: str):
    """
    Build nodes/edges for streamlit-agraph from ChemProp1 score dataframe.

    - Creates one node per unique cluster id in df (and ensures source/target are present)
    - Node label: m/z (2 decimals) when available, otherwise the cluster id
    - Node title: id, name, m/z, RT (best-effort)
    - Colors: source=blue, target=red, others=lightgray
    - Edges: direction based on Sign_ChemProp1 when abs_ChemProp1 > 0
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
    #    - direction based on Sign_ChemProp1 (1: a->b, -1: b->a)
    for _, r in df2.iterrows():
        a = r["CLUSTERID1"]
        b = r["CLUSTERID2"]

        abs_cp = r.get("abs_ChemProp1", 0)
        sign = r.get("Sign_ChemProp1", 0)

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


def get_inputs_from_state():
    required = ["nw", "chemprop1_ft", "chemprop1_md"]
    if not all(k in st.session_state for k in required):
        return None

    network_df = st.session_state["nw"].copy()
    features_df = st.session_state["chemprop1_ft"].copy()
    metadata_df = st.session_state["chemprop1_md"].copy()

    if network_df.empty or features_df.empty or metadata_df.empty:
        return None
    
    return network_df, features_df, metadata_df

def chemprop1_controls():
    if "run_chemprop1" not in st.session_state:
        st.session_state.run_chemprop1 = False

    show_options = st.checkbox("Run ChemProp1 score calculation", value=False)
    st.session_state.run_chemprop1 = bool(show_options)

    if not show_options:
        return False, None

    mode = st.radio(
        "Run ChemProp1 on:",
        ("Provided Edge Table", "Cascade edges", "User defined edge"),
        horizontal=True,
        key="chemprop1_mode"
    )

    return True, mode


def run_chemprop1_pipeline(network_df, features_df, metadata_df):
    
    chemprop1_df = run_chemprop1_analysis(network_df, features_df, metadata_df)

    if chemprop1_df is None or chemprop1_df.empty:
        return None

    chemprop1_df = drop_blank_score_rows(chemprop1_df, base_cols=network_df.shape[1])
    chemprop1_df = add_mz_rt_from_ft(chemprop1_df, st.session_state["ft"])

    # Add names if GNPS annotations exist
    if "an_gnps" in st.session_state and isinstance(st.session_state["an_gnps"], pd.DataFrame):
        gnps_df = st.session_state["an_gnps"]
        if not gnps_df.empty:
            chemprop1_df = add_names_from_gnps(chemprop1_df, gnps_df)

    st.session_state["ChemProp1_scores"] = chemprop1_df
    st.write(
        f"ChemProp1 Scoring Results " 
        f"({chemprop1_df.shape[0]} rows × {chemprop1_df.shape[1]} columns):"
        )
    st.dataframe(chemprop1_df, hide_index=True, use_container_width=True)
    return chemprop1_df


def drop_blank_score_rows(df, base_cols: int):
    # Keep everything in the original edge table + drop rows where all extra cols are NaN
    extra_cols = df.columns[base_cols:]
    if len(extra_cols) == 0:
        return df
    return df.dropna(subset=extra_cols, how="all")

##########################################
from collections import defaultdict, deque

def build_cascade_edges(
    edges: pd.DataFrame,
    source_id: int,
    max_neighbor: int = 15,
    col1: str = "CLUSTERID1",
    col2: str = "CLUSTERID2",
    component_col: str = "ComponentIndex",
) -> pd.DataFrame:
    """
    Build a cascade edge table from an undirected edge list.

    - Finds ComponentIndex group(s) containing `source_id`
    - Restricts to those components
    - BFS from source up to `max_neighbor` hops
    - Returns one row per reached node with hop distance in 'Neighbor'

    Output columns: CLUSTERID1 (source), CLUSTERID2 (node), Neighbor, ComponentIndex
    """

    required = {col1, col2, component_col}
    missing = required - set(edges.columns)
    if missing:
        raise ValueError(f"Edge table missing required columns: {sorted(missing)}")

    # Ensure numeric ids (but keep NaNs out)
    df = edges[[col1, col2, component_col]].copy()
    df = df.dropna(subset=[col1, col2, component_col])

    # Identify which components contain the source
    in_rows = df[(df[col1] == source_id) | (df[col2] == source_id)]
    if in_rows.empty:
        # Return empty with correct schema
        return pd.DataFrame(columns=[col1, col2, "Neighbor", component_col])

    source_components = sorted(in_rows[component_col].unique().tolist())

    # Restrict to those components
    df = df[df[component_col].isin(source_components)].copy()

    # Build adjacency per component
    # adj[component][node] -> set(neighbors)
    adj = defaultdict(lambda: defaultdict(set))
    for a, b, comp in zip(df[col1].tolist(), df[col2].tolist(), df[component_col].tolist()):
        adj[comp][a].add(b)
        adj[comp][b].add(a)

    # BFS per component, then combine results
    out_rows = []
    for comp in source_components:
        if source_id not in adj[comp]:
            continue

        visited = {source_id: 0}
        q = deque([source_id])

        while q:
            u = q.popleft()
            du = visited[u]
            if du == max_neighbor:
                continue
            for v in adj[comp].get(u, []):
                if v not in visited:
                    visited[v] = du + 1
                    q.append(v)

        # Add reached nodes except source itself
        for node, dist in visited.items():
            if node == source_id:
                continue
            if 1 <= dist <= max_neighbor:
                out_rows.append({
                    col1: source_id,          # source
                    col2: node,               # reached node
                    "Neighbor": dist,         # hop distance
                    component_col: comp
                })

    # Create final table (deduplicate in case a node is reachable in multiple components)
    out = pd.DataFrame(out_rows)
    if out.empty:
        return pd.DataFrame(columns=[col1, col2, "Neighbor", component_col])

    out = out.sort_values([component_col, "Neighbor", col2]).drop_duplicates(
        subset=[col1, col2, component_col], keep="first"
    ).reset_index(drop=True)

    return out

def render_results_summary(network_df, result_df):

    st.info(
        f"""
    **Edge Filtering Summary**

    - 📊 **Original edge table:** {network_df.shape[0]} edges  
    - 🧹 **ChemProp2 edge table:** {result_df.shape[0]} edges
    """
    )
  
    if network_df.shape[0] != result_df.shape[0]:
        st.warning("The reduced number of edges might be due to the removal of blank entries.")


def render_download_buttons(result_df):
    user_filename = st.text_input(
        "Enter the filename for the CSV and press Enter to apply the name:",
        value="chemprop1_scores_results.csv",
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
            generate_graphml_zip_chemprop1()

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


def render_scores_plot(scores_df):
    df = scores_df.copy()
    df["CLUSTER IDs"] = df["CLUSTERID1"].astype(str) + "-" + df["CLUSTERID2"].astype(str)

    fig = px.scatter(
        df,
        x="DeltaMZ",
        y="ChemProp1",
        hover_name="CLUSTER IDs",
        title="Scatter Plot: DeltaMZ vs ChemProp1 scores",
        labels={"DeltaMZ": "Delta M/Z", "ChemProp1": "ChemProp1 scores"},
    )
    st.plotly_chart(fig, use_container_width=True)

def render_filters_and_plots(scores_df, features_df, metadata_df):
    if scores_df is None or scores_df.empty:
        return

    if not st.checkbox("Show Filters", key="show_filters_chemprop"):
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

        filtered_df = get_filtered_edges_ui(df, filter_mode)
        st.dataframe(filtered_df, hide_index=True, use_container_width=True)
    render_edge_detail_plots(filtered_df, scores_df, features_df, metadata_df)
    if (st.session_state.get("gnps_task_id") or "").strip():
        render_spectra_modifinder(filtered_df, st.session_state.get("an_gnps"))


def get_filtered_edges_ui(df, filter_mode):
    c1, c2, c3 = st.columns(3)

    if filter_mode == "Filter by Score":
        with c1:
            score_min = st.number_input("Min Score:", value=float(df["ChemProp1"].min()))
        with c2:
            score_max = st.number_input("Max Score:", value=float(df["ChemProp1"].max()))

        sub = df[(df["ChemProp1"] >= score_min) & (df["ChemProp1"] <= score_max)].copy()

        with c3:
            if sub.empty:
                st.warning("No edges in this score range.")
                return sub
            sub["Dropdown_Display"] = (
                "ID: " + sub["CLUSTER IDs"].astype(str) + ", Score: " + sub["ChemProp1"].round(3).astype(str)
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


def render_edge_detail_plots(filtered_df, all_scores_df, features_df, metadata_df):
    if filtered_df is None or filtered_df.empty:
        return

    with st.container():
        c1, c2 = st.columns(2)

        with c1:
            selected_row = filtered_df.iloc[0]
            fig, status = plot_intensity_trends_single_row_chemprop1(
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

                allowed = ["ComponentIndex", "Cosine", "DeltaMZ", "ChemProp1", "Sign_ChemProp1", "abs_ChemProp1"]
                available = [c for c in plot_df.columns if c in allowed]

                edge_label_column = st.selectbox(
                    "Select column for edge labels:",
                    options=available,
                    help="To save the network image as PNG, right-click on empty space and select 'Save image as'.",
                )

                nodes, edges = generate_graph_from_df_chemprop1(plot_df, filtered_df, edge_label_column)

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


def _legend_html():
    return """
    <div style="display: flex; justify-content: center;">
        <div style="display: flex; align-items: center; margin-right: 20px;">
            <div style="width: 40px; height: 40px; border-radius: 50%; background-color: blue; color: white;
                        display: flex; align-items: center; justify-content: center; font-size: 12px;">
                Source
            </div>
            <div style="margin-top:10px; font-size:12px;">
                CLUSTERID1</b>
            </div>
        </div>
        <div style="display: flex; align-items: center;">
            <div style="width: 40px; height: 40px; border-radius: 50%; background-color: red; color: white;
                        display: flex; align-items: center; justify-content: center; font-size: 12px;">
                Target
            </div>
            <div style="margin-top:10px; font-size:12px;">
                CLUSTERID2</b>
            </div>
        </div>
    </div>
    """


def build_dashinterface_scan_url(task_id: str, scan1: str, scan2: str) -> str:
    usi1 = build_usi(task_id, scan1)
    usi2 = build_usi(task_id, scan2)

    params = {
        "usi1": usi1,
        "usi2": usi2,
        "width": 10.0,
        "height": 6.0,
        "cosine": "standard",
        "fragment_mz_tolerance": 0.1,
        "grid": True,
        "annotate_precision": 4,
        "annotation_rotation": 90,
    }
    return "https://metabolomics-usi.gnps2.org/dashinterface/?" + up.urlencode(
        params, safe=":/[](),"
    )


def build_usi(task_id: str, scan_id) -> str:
    return (
        f"mzspec:GNPS2:TASK-{task_id}-nf_output/clustering/"
        f"spectra_reformatted.mgf:scan:{scan_id}"
    )


def build_modifinder_url(
    usi1: str,
    usi2: str,
    *,
    smiles1: str = "",
    smiles2: str = "",
    adduct: str = "[M+H]1+",
    ppm_tolerance: int = 10,
    base_peak_filter_ratio: float = 0.01,
    helpers: str = "",
) -> str:
    base_url = "https://modifinder.gnps2.org/"

    params = {
        "USI1": usi1,
        "USI2": usi2,
        "SMILES1": smiles1 or "",
        "SMILES2": smiles2 or "",
        "Helpers": helpers,
        "adduct": adduct,
        "ppm_tolerance": ppm_tolerance,
        "filter_peaks_variable": base_peak_filter_ratio,
    }

    return f"{base_url}?{up.urlencode(params, quote_via=up.quote)}"



def get_smiles_from_annotation(an_gnps, scan_id):
    """
    Lookup SMILES for a given scan ID from GNPS annotation table.

    Parameters
    ----------
    an_gnps : pd.DataFrame or None
        GNPS annotation table (must contain '#Scan#' and 'smiles' columns)
    scan_id : str or int

    Returns
    -------
    str
        SMILES string if found, else empty string
    """
    if an_gnps is None or an_gnps.empty:
        return ""

    if "#Scan#" not in an_gnps.columns or "Smiles" not in an_gnps.columns:
        return ""

    scan_id = str(scan_id)

    hit = an_gnps.loc[
        an_gnps["#Scan#"].astype(str) == scan_id, "Smiles"
    ]

    if hit.empty:
        return ""

    smiles = hit.iloc[0]
    return "" if pd.isna(smiles) else str(smiles).strip()

def render_spectra_modifinder(filtered_df, an_gnps, task_id_key: str = "gnps_task_id"):
    if filtered_df is None or filtered_df.empty:
        return

    task_id = (st.session_state.get(task_id_key) or "").strip()
    if not task_id:
        return  # hard guard (no UI)

    row = filtered_df.iloc[0]
    scan1 = str(row.get("CLUSTERID1", "")).strip()
    scan2 = str(row.get("CLUSTERID2", "")).strip()

    if not (scan1 and scan2):
        return

    # --- Build USIs
    usi1 = build_usi(task_id, scan1)
    usi2 = build_usi(task_id, scan2)

    # --- Lookup SMILES (optional)
    smiles1 = get_smiles_from_annotation(an_gnps, scan1)
    smiles2 = get_smiles_from_annotation(an_gnps, scan2)

    dash_url = build_dashinterface_scan_url(task_id, scan1, scan2)

    modifinder_url = build_modifinder_url(
        usi1=usi1,
        usi2=usi2,
        smiles1=smiles1,
        smiles2=smiles2,
    )

    c1, c2 = st.columns(2)
    with c1:
        st.link_button(
            f"🔎 Spectrum Resolver (scan {scan1} vs {scan2})",
            dash_url,
            type="primary",
            use_container_width=True,
        )

    with c2:
        st.link_button(
            "🧩 Open in ModiFinder",
            modifinder_url,
            use_container_width=True,
        )

        if not smiles1 and not smiles2:
            st.caption("No SMILES found in GNPS annotations (can be added manually in ModiFinder).")

