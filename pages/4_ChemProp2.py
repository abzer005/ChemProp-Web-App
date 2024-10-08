# Import necessary libraries
import streamlit as st
from src.chemprop1_helpers import *
import pandas as pd
import numpy as np
import re

st.markdown("### Under construction")
st.markdown("## ChemProp2 Analysis")

st.markdown("""
            **Instructions:** 
            For ChemProp2 analysis, select an attribute from your dataset that contains more than two distinct sequential points (e.g., multiple time points or spatial points). 
            This analysis focuses on the correlation relationship between pairs of features. It operates on the principle that if two features, 
            say A and B, show a positive correlation (both increasing or decreasing together over time or space), they will receive a score close to zero. 
            Conversely, if they exhibit anti-correlation (A increases while B decreases, or vice versa), they will receive a score reflecting this relationship. 
            The resulting ChemProp2 score ranges from -1 to +1, with the magnitude indicating the strength of the potential biotransformation and the sign suggesting the direction 
            (a positive score implies transformation from A to B, while a negative score implies transformation from B to A).
            """)

# Check if metadata is loaded
if 'md_for_analysis' in st.session_state and not st.session_state['md_for_analysis'].empty:
    
    chemprop_md = st.session_state['md_for_analysis']

    df = inside_levels(chemprop_md)
    mask = df.apply(lambda row: len(row['LEVELS']) == 0, axis=1)
    df = df[~mask]
    
    # Display the overview dataframe and checkboxes
    st.markdown("### Select Metadata Columns for ChemProp2 Calculation")

    # Add a new column for checkboxes (initialized to False)
    df['ChemProp2'] = False
    df['Treatments'] = False

    # Display an editable DataFrame using st.data_editor
    edited_df = st.data_editor(df, num_rows="fixed", 
                               use_container_width=True, 
                               hide_index = True,
                               disabled=("ATTRIBUTES", "LEVELS", "COUNTS"),
                               column_order=("ChemProp2", "Treatments" , "ATTRIBUTES", "LEVELS", "COUNTS")
                               )

    # Collect the rows where checkboxes are selected (True)
    chemprop2_row = edited_df[edited_df['ChemProp2'] == True]['ATTRIBUTES'].tolist()
    treatment_row = edited_df[edited_df['Treatments'] == True]['ATTRIBUTES'].tolist()

    # Show the selected attributes for Treatments
    if chemprop2_row:
        st.success(f"The column selected for calculating ChemProp2: {', '.join(chemprop2_row)}")
    else:
        st.warning("No column selected for calculating ChemProp2 scores.")


    # Show the selected attributes for Treatments
    if treatment_row:
        st.success(f"The column selected for Treatment conditions: {', '.join(treatment_row)}")
    else:
        st.warning("No column selected for Treatments.")
    
    # Subset the DataFrames based on the user's selection
    treatment_md = chemprop_md[treatment_row] if treatment_row else pd.DataFrame()
    time_md = chemprop_md[chemprop2_row]

    # Apply the function to the subset DataFrame
    time_md = time_md.applymap(strip_non_numeric)

    # Convert the stripped values to numeric types (float)
    time_md = time_md.apply(pd.to_numeric, errors='coerce')

    # Conditional creation of subset_chemprop_md
    if treatment_row:
        # If both time and treatment columns are selected, concatenate them
        subset_chemprop_md = pd.concat([time_md, treatment_md], axis=1, join='outer')
    else:
        # If only time columns are selected, use time_md
        subset_chemprop_md = time_md
        
    st.markdown("### ChemProp2 Metadata with Selected Columns:")
    st.dataframe(subset_chemprop_md)

    final_ft = st.session_state['ft_for_analysis'].copy()
    st.dataframe(final_ft)
    
    st.session_state['chemprop_ft'] = final_ft
    st.session_state['chemprop_subset_md'] = subset_chemprop_md
    st.session_state['chemprop_md'] = time_md

    subset_md = st.checkbox("Do you want to subset the data?", False)
    if subset_md:
        st.markdown("Select samples based on their treatment conditions.")
        # Sample and blank selection interface
        c1, c2 = st.columns(2)
        
        # Let the user select the attribute column for categorization
        attribute_column = c1.selectbox("Select the column for categorization", subset_chemprop_md.columns)

        # Allow the user to select multiple categories as a group
        subset_group = c2.multiselect("Select categories for filtering", options=sorted(subset_chemprop_md[attribute_column].astype(str).unique()))
        subset_ft = st.session_state['ft_for_analysis'].copy()

        if subset_group:
            subgroup_indices = subset_chemprop_md[subset_chemprop_md[attribute_column].isin(subset_group)].index
            subgroup = subset_ft.loc[:, subgroup_indices]
            subgroup_md = subset_chemprop_md.loc[subgroup_indices]

            with st.expander(f"Selected Group {subgroup.shape}"):
                st.dataframe(subgroup)
                st.dataframe(subgroup_md)
                st.session_state['chemprop_ft'] = subgroup
                st.session_state['chemprop_subset_md'] = subgroup_md
                st.session_state['chemprop_md'] = subgroup_md[chemprop2_row]
        else:
            st.warning("No group selected.")
         
else:
    st.warning("Metadata not loaded in the session state.")



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

