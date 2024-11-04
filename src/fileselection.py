import streamlit as st
from .common import *  # Importing common functionalities from the 'common' module
import pandas as pd
from gnpsdata import taskresult
from gnpsdata import workflow_fbmn
import urllib

patterns = [
    ["m/z", "mz", "mass over charge"],
    ["rt", "retention time", "retention-time", "retention_time"],
]

allowed_formats = "Allowed formats: csv (comma separated), tsv (tab separated), txt (tab separated), xlsx (Excel file)."

def string_overlap(string, options):
    """
    Check if any of the given options are present in the string.
    Exclude any string containing "mzml".

    Parameters:
    string (str): The string to be checked.
    options (list): A list of substrings to search for in the string.

    Returns:
    bool: True if any option is found in the string and "mzml" is not present, False otherwise.
    """
    for option in options:
        if option in string and "mzml" not in string:
            return True
    return False

def load_example():
    """
    Load example datasets into Streamlit's session state.
    """
    # Reset session state data
    for key in ['ft', 'md', 'nw', 'an_gnps', 'an_analog']:
        st.session_state[key] = None
        
    st.session_state['ft'] = open_df("example-data/FeatureMatrix.csv").set_index("row ID")
    st.session_state['md'] = open_df("example-data/MetaData.txt").set_index("filename")
    st.session_state['nw'] = open_df("example-data/NetworkNodePairs.csv")
    st.session_state['an_gnps'] = open_df("example-data/GNPSannotations.tsv")

    #st.session_state['ft'] = open_df("example-data/dummy_featuretable.csv").set_index("row ID")
    #st.session_state['md'] = open_df("example-data/dummy_metadata.txt").set_index("filename")
    #st.session_state['nw'] = open_df("example-data/dummy_edgetable.csv")
    #st.session_state['an_gnps'] = open_df("example-data/GNPSannotations.tsv")


#@st.cache_data  # Corrected cache decorator
#def get_networkpairs_dataframe(task, gnps2=True):
#    if gnps2:
#        return taskresult.get_gnps2_task_resultfile_dataframe(task, "nf_output/networking/filtered_pairs.tsv")
#    else:
#        return taskresult.get_task_resultview_dataframe(task, "filtered_pairs/", output_file)

#@st.cache_data  # Corrected cache decorator
#def get_gnpsannotations_dataframe(task, gnps2=True):
#    if gnps2:
#        return taskresult.get_gnps2_task_resultfile_dataframe(task, "nf_output/library/merged_results_with_gnps.tsv")
#    else:
#        return taskresult.get_task_resultview_dataframe(task, "merged_results_with_gnps/", output_file)


def load_ft(ft_file):
    """
    Load and process the feature table.

    Parameters:
    ft_file (file): The feature table file.

    Returns:
    DataFrame: Processed feature table.
    """
    ft = open_df(ft_file)
    ft = ft.dropna(axis=1)  # Drop columns with missing values
    return ft

def load_md(md_file):
    """
    Load and process metadata. Set 'filename' as the index if present.

    Parameters:
    md_file (file): The metadata file.

    Returns:
    DataFrame: Processed metadata.
    """
    md = open_df(md_file)
    return md

def load_nw(network_file):
    """
    Load and process network pair file. 
    """
    nw = open_df(network_file)
    return nw

def load_annotation(annotation_file):
    """
    Load and process annotation file 
    """
    an_gnps = open_df(annotation_file)
    return an_gnps

def display_dataframe_with_toggle(df_key, display_name):
    if df_key in st.session_state and isinstance(st.session_state[df_key], pd.DataFrame):
        st.write(f"### {display_name}")

        col1, col2 = st.columns([0.8, 0.2])

        # Show dimensions
        num_rows, num_cols = st.session_state[df_key].shape
        col1.write(f"Dimension: {num_rows} rows × {num_cols} columns")

        view_all = col2.checkbox("View all", key=f"{df_key}_toggle")

        if view_all:
            st.dataframe(st.session_state[df_key])  # Show full dataframe
        else:
            st.dataframe(st.session_state[df_key].head())  # Show header


def load_from_gnps(task_id):

    try: # GNPS2 will run here
        ft = workflow_fbmn.get_quantification_dataframe(task_id, gnps2=True).set_index("row ID")
        md = workflow_fbmn.get_metadata_dataframe(task_id, gnps2=True).set_index("filename")
        an = taskresult.get_gnps2_task_resultfile_dataframe(task_id, "nf_output/library/merged_results_with_gnps.tsv")
        nw = taskresult.get_gnps2_task_resultfile_dataframe(task_id, "nf_output/networking/filtered_pairs.tsv")

    except urllib.error.HTTPError as e:
        print(f"HTTP Error encountered: {e}") # GNPS1 task IDs can not be retrieved and throw HTTP Error 500
        
    if md.empty: # Handle empty metadata
        md = pd.DataFrame()

    return ft, md, an, nw
