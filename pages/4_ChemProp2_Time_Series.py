# Import necessary libraries
import streamlit as st
from src.chemprop1_helpers import *
import pandas as pd
import numpy as np
import re
import plotly.express as px


st.markdown("### Under construction")
st.markdown("## ChemProp2 Analysis for Time Series Data")

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

    st.markdown("""
                - ChemProp2: Select the metadata column that contains the timepoint information. This column must be present in the metadata for ChemProp2 calculation.
                - Treatments: Select the metadata column that contains the treatment conditions (e.g., control and treatment groups). This category is optional.
                """)

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
        
    st.markdown("### ChemProp2 Metadata with Selected Columns")
    st.dataframe(subset_chemprop_md)

    final_ft = st.session_state['ft_for_analysis'].copy()
    st.dataframe(final_ft)
    
    st.session_state['chemprop_ft'] = final_ft
    st.session_state['chemprop_subset_md'] = subset_chemprop_md
    st.session_state['chemprop_md'] = time_md
else:
    st.warning("Metadata not loaded in the session state.")


st.markdown("### Filter the data")
subset_md = st.checkbox("Do you want to subset the data further?", False)
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
         



# Streamlit Interface
st.markdown('## ChemProp2 Scores')

if 'nw' in st.session_state and 'chemprop_ft' in st.session_state and 'chemprop_md' in st.session_state:

    network_df = st.session_state['nw'].copy()
    features_df = st.session_state['chemprop_ft'].copy()
    metadata_df = st.session_state['chemprop_md'].copy()

    # Check if DataFrames are non-empty before processing
    if not network_df.empty and not features_df.empty and not metadata_df.empty:

        # Add a button to trigger ChemProp2 Scoring
        if st.button("Run ChemProp2"):
            node_names = list(features_df.index)
    
            result_df = ChemProp2_scoring(network_df, features_df, node_names, metadata_df)
        
            st.write("ChemProp2 Scoring Results:")

            # Get the number of columns in network_df
            network_cols = network_df.shape[1]

            # Define the range of columns in result_df you want to check (from 6th column to 11th column)
            extra_cols_range = result_df.columns[network_cols:]  # From ncol(network_df) + 1 to ncol(result_df)

            # Drop rows where all values in the selected columns are NaN
            result_df = result_df.dropna(subset=extra_cols_range, how='all')
            st.dataframe(result_df)
            st.session_state['ChemProp2_scores'] = result_df

            # Allow users to download the results
            # Let the user input the filename (default value provided)
            user_filename = st.text_input("Enter the filename for the CSV:", value="chemprop2_scores_results.csv")

            # Ensure the file ends with '.csv'
            if not user_filename.endswith('.csv'):
                user_filename += '.csv'
                
            # Convert result_df to CSV and encode it
            chemprop2_result = result_df.to_csv(index=False).encode('utf-8')
    
            # Provide a download button using the user-provided filename
            st.download_button(
                label="Download Results as CSV",
                data=chemprop2_result,
                file_name=user_filename,
                mime='text/csv',
                )
            
            #Getting the scores plot:
            scores_df = st.session_state['ChemProp2_scores'].copy()
            scores_df['CLUSTER IDs'] = scores_df['CLUSTERID1'].astype(str) + '-' + scores_df['CLUSTERID2'].astype(str)
            
            # User selects which column to use as the y-axis
            y_axis_option = st.selectbox(
                'Choose the column of scores',
                options=['ChemProp2', 'ChemProp_spearman', 'ChemProp_log', 'ChemProp_sqrt']
                )

            # Plot a scatter plot using Plotly
            fig = px.scatter(
                scores_df, 
                x='DeltaMZ', 
                y=y_axis_option, 
                hover_name='CLUSTER IDs',  # Show the 'Labels' column when hovering over a point
                title=f"Scatter Plot: DeltaMZ vs {y_axis_option}",
                labels={'DeltaMZ': 'Delta M/Z', y_axis_option: y_axis_option},
                range_y=[-1, 1]  # Set the y-axis range to be between -1 and +1
                )

            # Show the plot in Streamlit
            st.plotly_chart(fig)


        