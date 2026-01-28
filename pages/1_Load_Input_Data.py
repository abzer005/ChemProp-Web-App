# Import necessary libraries
import streamlit as st
import pandas as pd
from src.common import *        # Importing common functionalities
from src.fileselection import * # Importing file selection functionalities
from src.cleanup import *

page_setup()

# Introduction Section
st.markdown("### Please select your method for data input below.")

input_method = st.selectbox(
    "Select Input Method",
    ["Use Example Dataset", "Manual Input", "GNPS(2) FBMN task ID"],
    index=0,
    key="input_method",
)

# --- Initialize session_state keys ---
for key in ["ft", "md", "nw", "an_gnps"]:
    if key not in st.session_state:
        st.session_state[key] = None

if "gnps_task_id" not in st.session_state:
    st.session_state["gnps_task_id"] = ""

if "last_input_method" not in st.session_state:
    st.session_state["last_input_method"] = input_method

# --- Detect change in input method ---
# if st.session_state["last_input_method"] != input_method:
#     for key in ["ft", "md", "nw", "an_gnps"]:
#         st.session_state[key] = None

#     if st.session_state["last_input_method"] == "GNPS(2) FBMN task ID":
#         st.session_state["gnps_task_id"] = ""
        
#     st.session_state["last_input_method"] = input_method
#     st.rerun()


if st.session_state["last_input_method"] != input_method:
    KEEP_KEYS = {
        "last_input_method",
        "input_method",
        "page",
    }

    for key in list(st.session_state.keys()):
        if key not in KEEP_KEYS:
            del st.session_state[key]

    st.session_state["last_input_method"] = input_method
    st.rerun()

# ============================================================
# 1. Example Dataset (special: loads all 4 tables)
# ============================================================

if input_method == "Use Example Dataset":
    load_example()
    show_input_tables_in_tabs()
    #st.stop()  # don't show uploaders below in this mode

# ============================================================
# 2. Manual Input
# ============================================================
elif input_method == "Manual Input":
    st.info("💡 Upload tables in txt (tab separated), tsv, csv or xlsx (Excel) format.")

    # Create 2 columns for the ft, md file uploaders
    col1, col2 = st.columns(2)
    with col1:
        ft_file = st.file_uploader("Upload Feature Table", 
                                   type=["csv", "xlsx", "txt", "tsv"],
                                   help = "This table is a key output of LC-MS/MS metabolomics studies. The table presents a list of mass spectral features along with their relative intensities (represented by its integrated peak area) observed across various samples.")
        if ft_file:
            st.session_state['ft'] = load_ft(ft_file).set_index("row ID")

    with col2:
        md_file = st.file_uploader("Upload Metadata", 
                                   type=["csv", "xlsx", "txt", "tsv"],
                                   help = "The metadata table is created by the user, providing additional context for the measured samples, such as sample type, species, and tissue type, etc.")
        if md_file:
            st.session_state['md'] = load_md(md_file).set_index("filename")
    
    # Create 2 columns for the nw, annotation file uploaders
    col3, col4 = st.columns(2)
    with col3:
        network_file = st.file_uploader("Upload Edge File from FBMN", 
                                        type=["csv", "xlsx", "txt", "tsv"],
                                        help = "This file maps connections between node pairs in a metabolomic network. Each row represents a node pair, identified by Cluster IDs 1 and 2, corresponding to features in the quantification table. The file includes parameters like cosine score and mz difference, which could be used as edge connecting the nodes.")
        if network_file:
            st.session_state['nw'] = load_nw(network_file)
    
    with col4:
        annotation_file = st.file_uploader("Upload Annotation Information File (Optional)", type=["csv", "xlsx", "txt", "tsv"])
        if annotation_file:
            st.session_state['an_gnps'] = load_annotation(annotation_file)

    show_input_tables_in_tabs()

# ============================================================
# 3. GNPS(2) FBMN task ID
# ============================================================
elif input_method == "GNPS(2) FBMN task ID":
    
    st.warning("💡 This tool only supports task IDs from GNPS2.")

    task_id_default = "1b5b94b4191d4223a5f57afb2aaaf0b0"
    task_id = st.text_input("GNPS task ID", task_id_default)

    st.session_state["gnps_task_id"] = task_id.strip()

    _, c2, _ = st.columns(3)
    load_button_clicked = c2.button(
        "Load files from GNPS",
        type="primary",
        disabled=len(task_id) == 0,
        use_container_width=True
    )
    #give a default task ID for testing: 1b5b94b4191d4223a5f57afb2aaaf0b0
  
    if load_button_clicked and task_id:
        st.session_state["ft"], st.session_state["md"], st.session_state["an_gnps"], st.session_state["nw"] = load_from_gnps(task_id)

        if st.session_state.get("ft") is not None and st.session_state.get("nw") is not None and (
            st.session_state.get("md") is None or st.session_state["md"].empty
        ):
            st.warning("Metadata is missing. Please upload it manually.")
            md_file = st.file_uploader("Upload Metadata", type=["csv", "xlsx", "txt", "tsv"])
            if md_file:
                st.session_state["md"] = load_md(md_file)
                st.success("Metadata was loaded successfully!")

    show_input_tables_in_tabs()


st.markdown("## Data Cleanup")

# Check if the data is available in the session state
if (
    'ft' in st.session_state and 
    'md' in st.session_state and 
    st.session_state['ft'] is not None and 
    not st.session_state['ft'].empty and 
    st.session_state['md'] is not None and 
    not st.session_state['md'].empty
):

    ft = st.session_state['ft'].copy()
    md = st.session_state['md'].copy()

    # Initialize empty DataFrame for samples metadata
    samples_md = pd.DataFrame()
    samples = pd.DataFrame()

    # If data is available, proceed with cleanup and checks
    cleaned_ft = clean_up_ft(ft)
    cleaned_md = clean_up_md(md)

    # Check if ft column names and md row names are the same
    cleaned_md, cleaned_ft = check_columns(cleaned_md, cleaned_ft)
    st.markdown("#### Metadata overview")
    df = inside_levels(cleaned_md)
    mask = df.apply(lambda row: len(row['LEVELS']) == 0, axis=1)
    df = df[~mask]
    st.dataframe(df)

    st.session_state['ft_for_analysis'] = cleaned_ft
    st.session_state['md_for_analysis'] = cleaned_md

    blank_removal = st.checkbox("Remove blank features?", False)
    if blank_removal:
         # Select true sample files (excluding blank and pools)
         st.markdown("#### Samples Selection")
         st.markdown("Select samples (excluding blank and pools) based on the above table.")
         # Sample and blank selection interface
         c1, c2, c3 = st.columns(3)
         
         # Let the user select the attribute column for categorization
         attribute_column = c1.selectbox("Select the column for categorization", cleaned_md.columns)
         
         # Allow the user to select multiple categories as samples
         sample_categories = c2.multiselect("Select categories for samples", options=sorted(cleaned_md[attribute_column].astype(str).unique()))

         # Allow the user to select multiple categories as blanks
         blank_categories = c3.multiselect("Select categories for blanks", options=sorted(cleaned_md[attribute_column].astype(str).unique()))

         # Filter the dataframes based on the selections
         if sample_categories:
            sample_indices = cleaned_md[cleaned_md[attribute_column].isin(sample_categories)].index
            samples = cleaned_ft.loc[:, sample_indices]
            samples_md = cleaned_md.loc[sample_indices]

            with st.expander(f"Selected samples {samples.shape}"):
                st.dataframe(samples)
                st.session_state['samples'] = samples

         else:
            st.warning("No samples selected.")
    
         if blank_categories:
            blank_indices = cleaned_md[cleaned_md[attribute_column].isin(blank_categories)].index
            blanks = cleaned_ft.loc[:, blank_indices]

            with st.expander(f"Selected blanks {blanks.shape}"):
                st.dataframe(blanks)
                st.session_state['blanks'] = blanks
            
         else:
            st.warning("No blanks selected.")
            
         if samples.shape[1] != ft.shape[1] and blank_categories:
             c1, c2 = st.columns(2)
             cutoff = c1.number_input(
                 "cutoff threshold for blank removal",
                 min_value=0.01,
                 max_value=1.0,
                 value=0.3,
                 step=0.05,
                 help="""The recommended cutoff range is between 0.1 and 0.3.
                 Features with intensity ratio of (blank mean)/(sample mean) above the threshold (e.g. 30%) are considered noise/background features.
                 """)
             ft, n_background_features, n_real_features = remove_blank_features(blanks, samples, cutoff)
             c2.metric("background or noise features", n_background_features)

             with st.expander(f"Feature table after removing blanks {ft.shape}"):
                 show_table(ft, "blank-features-removed")
                 
             # Update the feature table in the session state after removing blanks
             st.session_state['ft_for_analysis'] = ft
             st.session_state['md_for_analysis'] = samples_md

         else:
             st.warning("You selected everything as sample type. Blank removal is not possible.")
         
         st.session_state['blank_removal_done'] = True
    else:
         st.session_state['blank_removal_done'] = False    
    
    
    if not st.session_state['ft_for_analysis'].empty:
        cutoff_LOD = get_cutoff_LOD(st.session_state['ft_for_analysis'])

    imputation = st.checkbox("Impute missing values?", 
                             False, 
                             help=f"These values will be filled with random number between 0.1 and {cutoff_LOD} (Limit of Detection) during imputation.")
    if imputation:
        if cutoff_LOD > 0.1:
            imputed_ft = impute_missing_values(st.session_state['ft_for_analysis'], cutoff_LOD)
            with st.expander(f"Imputed data - features: {imputed_ft.shape[0]}, samples: {imputed_ft.shape[1]}"):
                st.dataframe(imputed_ft)
            
            st.session_state['ft_for_analysis'] = imputed_ft
        else:
            st.warning(f"Can't impute with random values between 1 and lowest value, which is {cutoff_LOD} (rounded).")
        
        st.session_state['imputation_done'] = True
    else:
        st.session_state['imputation_done'] = False


    perform_normalization = st.checkbox("Perform TIC normalization?", False)
    if perform_normalization:
        normalized_md, normalized_ft = normalization(st.session_state['ft_for_analysis'], 
                                                     st.session_state['md_for_analysis'])
        normalized_ft = normalized_ft.T
        
        with st.expander(f"Normalized data - features: {normalized_ft.shape[0]}, samples: {normalized_ft.shape[1]}"):
            st.dataframe(normalized_ft)
            
        st.session_state['ft_for_analysis'] = normalized_ft

        st.session_state['normalization_done'] = True
    else:
        st.session_state['normalization_done'] = False
else:
    # If data is not available, display a message
    st.warning("Data not loaded. Please load the data first.")

st.warning(
    "After data cleanup, proceed to the **ChemProp1** or **ChemProp2** analysis pages "
    "based on the number of timepoints in your dataset. "
    "Use **ChemProp1** for datasets with **exactly two timepoints**, and **ChemProp2** "
    "for datasets with **three or more timepoints**. "
    "You can access both analysis modules from the **left-hand sidebar navigation menu**, "
    "where each module is listed as a separate page."
)

# Displaying the message with a link
st.markdown("""
    <style>
    .big-font {
        font-size:20px !important;
    }
    </style>
    <div class='big-font'>
    For more functionalities such as univariate and multivariate statistics, 
    explore our <a href="https://fbmn-statsguide.gnps2.org/" target="_blank"><b>FBMN-STATS webpage</b></a>.
    </div>
    """, unsafe_allow_html=True)
