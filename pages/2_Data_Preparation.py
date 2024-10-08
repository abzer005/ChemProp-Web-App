#Import necessary libraries
import streamlit as st
from src.common import *
from src.fileselection import *
from src.cleanup import *
import pandas as pd
import numpy as np

st.markdown("## Blank removal")

# Check if the data is available in the session state
if 'ft' in st.session_state and 'md' in st.session_state:
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
                 0.1, 1.0, 0.3, 0.05,
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

else:
    # If data is not available, display a message
    st.warning("Data not loaded. Please go back and load the data first.")


# Displaying the message with a link
st.markdown("""
    <style>
    .big-font {
        font-size:20px !important;
    }
    </style>
    <div class='big-font'>
    For more functionalities such as data cleanup, univariate and multivariate statistics, 
    explore our <a href="https://fbmn-statsguide.gnps2.org/" target="_blank"><b>FBMN-STATS webpage</b></a>.
    </div>
    """, unsafe_allow_html=True)

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
            df_dimensions.append(f"{df.shape}")

    # Display the table if there are dataframes
    if df_names:
        # Convert lists to pandas DataFrame for display
        df_table = pd.DataFrame({
            'Dataframe': df_names,
            'Dimensions (rows x cols)': df_dimensions
        })
        st.table(df_table)