# Import necessary libraries
import streamlit as st
from src.chemprop1_helpers import *
import pandas as pd
import numpy as np

st.markdown("## Under construction")
st.markdown("### ChemProp1 Analysis")

st.markdown("""
**Instructions:** 
For ChemProp1 analysis, select an attribute from your dataset that contains exactly two distinct sequential points (e.g., two different time points or two spatial points). 
This analysis is based on the assumption that the abundance of an initial compound decreases over time or across space, while the abundance of a new compound increases.
""")

# Check if metadata is loaded
if 'md_for_analysis' in st.session_state and not st.session_state['md_for_analysis'].empty:
    st.markdown("#### Metadata Table:")
    df = inside_levels(st.session_state.md_for_analysis)
    mask = df.apply(lambda row: len(row['LEVELS']) == 0, axis=1)
    df = df[~mask]
    st.dataframe(df)

    # Generate the list of attributes that contain numerical values
    chemprop_md = st.session_state.md_for_analysis.copy()
    valid_attributes = [c for c in chemprop_md.columns if contains_convertible_values(chemprop_md[c])]
    c1, c2 = st.columns(2)

    c1.radio("Select the attribute that is appropriate for ChemProp1 Calculation:",
        valid_attributes,
        key = "chemprop_attr",
       horizontal=False,)
    st.write("You selected:", st.session_state.chemprop_attr)
    unique_values = chemprop_md[st.session_state.chemprop_attr].dropna().unique()

    if len(unique_values) < 3:
        st.success(f"You have only {len(unique_values)} timepoints, hence using ChemProp1 is appropriate for this data.", icon="✅")
        #st.write(f"The values are: {', '.join(str(value) for value in unique_values)}")
    else:
        st.error("You cannot run ChemProp1 for timepoints more than 2. Please run ChemProp2.", icon="🚨")

    # Column 2: Toggle to remove characters from the metadata
    on = c2.toggle('Strip away characters from the attribute')
    if on:
        chemprop_md[st.session_state.chemprop_attr] = convert_series_values(chemprop_md[st.session_state.chemprop_attr])
        st.session_state['chemprop_md'] = chemprop_md
        c2.write('Characters removed!')
        converted_unique_values = chemprop_md[st.session_state.chemprop_attr].unique()
        c2.write(f"The new values are: {', '.join(str(value) for value in converted_unique_values)}")
    else:
        # Reset to original data if toggle is off
        st.session_state['chemprop_md'] = st.session_state.md_for_analysis.copy()
        c2.write("Reset to original data.")
else:
    st.warning("Metadata not loaded. Please load the metadata first.")

subset_md = st.checkbox("Subset Metadata?", False)

# Check if the required data is available in the session state
if 'ft_for_analysis' in st.session_state and 'chemprop_md' in st.session_state and 'nw' in st.session_state:
    # Determine which data to use for ChemProp1 analysis
    if 'subset_md_after' in st.session_state and 'subset_ft_after' in st.session_state and not st.session_state['subset_md_after'].empty and not st.session_state['subset_ft_after'].empty:
        md_to_use = st.session_state['subset_md_after'].copy()
        ft_to_use = st.session_state['subset_ft_after'].copy()
    else:
        md_to_use = st.session_state['chemprop_md'].copy()
        ft_to_use = st.session_state['ft_for_analysis'].copy()

    # Ask the user for confirmation to proceed with ChemProp1 analysis
    if st.button('Run ChemProp1'):
        st.session_state['nw_with_scores'] = run_chemprop1_analysis(st.session_state['nw'], 
                                               md_to_use, 
                                               ft_to_use, 
                                               st.session_state['chemprop_attr'])

        # Display results in tabs
        t1, t2 = st.tabs(["Network Table with Scores", "Scatter Plot with Scores"])
        with t1:
            st.dataframe(st.session_state['nw_with_scores'])
        with t2:
            st.write("Scatter plot functionality to be implemented")
            # Code for creating and displaying scatter plot goes here
else:
    st.warning("Required data for analysis is not loaded. Please check your dataset.")
