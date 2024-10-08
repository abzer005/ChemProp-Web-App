# Import necessary libraries
import streamlit as st
import pandas as pd
from src.common import *        # Importing common functionalities
from src.fileselection import * # Importing file selection functionalities

# Introduction Section
st.markdown("### Please select your method for data input below.")

# Input Selection Section
input_method = st.selectbox("Select Input Method", 
                            ["Use Example Dataset",
                             "Manual Input"
                             ])

# Clearing the session state 
if 'last_input_method' not in st.session_state:
    st.session_state['last_input_method'] = None

# Example Dataset Section
elif input_method == "Use Example Dataset":
    # Check if input method has changed
    if st.session_state['last_input_method'] != input_method:
        # Clear the data
        for key in ['ft', 'md', 'nw', 'an_gnps']:
            st.session_state[key] = None

        # Update the last input method
        st.session_state['last_input_method'] = input_method
    
    load_example()  # Load data into session state

    for file_name, key in zip(["Feature Matrix", "MetaData", "Network Node Pairs", "GNPS Annotations"],
                              ['ft', 'md', 'nw', 'an_gnps']):
        display_dataframe_with_toggle(key, file_name)

# Manual Input Section
if input_method == "Manual Input":
    if st.session_state['last_input_method'] != input_method:
        # Clear the data
        for key in ['ft', 'md', 'nw', 'an_gnps']:
            st.session_state[key] = None
        # Update the last input method
        st.session_state['last_input_method'] = input_method

    st.info("💡 Upload tables in txt (tab separated), tsv, csv or xlsx (Excel) format.")

    # Create 2 columns for the ft, md file uploaders
    col1, col2 = st.columns(2)
    with col1:
        ft_file = st.file_uploader("Upload Feature Table", 
                                    type=["csv", "xlsx", "txt", "tsv"],
                                    help = "This table is a key output of LC-MS/MS metabolomics studies. The table presents a list of mass spectral features along with their relative intensities (represented by its integrated peak area) observed across various samples.")
        if ft_file:
            st.session_state['ft'] = load_ft(ft_file)

    with col2:
        md_file = st.file_uploader("Upload Metadata", 
                                   type=["csv", "xlsx", "txt", "tsv"],
                                   help = "The metadata table is created by the user, providing additional context for the measured samples, such as sample type, species, and tissue type, etc.")
        if md_file:
            st.session_state['md'] = load_md(md_file)
    
    # Create 2 columns for the nw, annotation file uploaders
    col3, col4 = st.columns(2)
    with col3:
        network_file = st.file_uploader("Upload Network Edge File from FBMN", 
                                        type=["csv", "xlsx", "txt", "tsv"],
                                        help = "This file maps connections between node pairs in a metabolomic network. Each row represents a node pair, identified by Cluster IDs 1 and 2, corresponding to features in the quantification table. The file includes parameters like cosine score and mz difference, which could be used as edge connecting the nodes.")
        if network_file:
            st.session_state['nw'] = load_nw(network_file)
    
    with col4:
        annotation_file = st.file_uploader("Upload Annotation Information File (Optional)", type=["csv", "xlsx", "txt", "tsv"])
        if annotation_file:
            st.session_state['an_gnps'] = load_annotation(annotation_file)

    # Display headers and 'View all' buttons for each file
    for key, label in zip(['ft', 'md', 'nw', 'an_gnps'],
                          ["Feature Table", "Metadata", "Network Edge File", "Annotation Information File"]):
        
        if key in st.session_state and st.session_state[key] is not None:
            df = st.session_state[key]
            col1, col2 = st.columns([0.8, 0.2])
            col1.write(f"{label} - Header")
            view_all = col2.checkbox("View all", key=f"{key}_toggle")
            
            if view_all:
                st.dataframe(df)
            else:
                st.dataframe(df.head())

# Displaying Dataframes in Sidebar
with st.sidebar:
    st.write("## Uploaded Data Overview")

    # Create lists for dataframe information
    df_names = []
    df_dimensions = []

    # Iterate over session state items
    for df_name in st.session_state.keys():
        # Check if the item is a DataFrame
        if isinstance(st.session_state[df_name], pd.DataFrame):
            df = st.session_state[df_name]
            df_names.append(df_name)
            df_dimensions.append(f"{df.shape[0]} rows × {df.shape[1]} columns")

    # Display the table if there are dataframes
    if df_names:
        # Convert lists to pandas DataFrame for display
        df_table = pd.DataFrame({
            'Dataframe': df_names,
            'Dimensions': df_dimensions
        })
        st.table(df_table)