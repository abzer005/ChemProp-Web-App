# Import necessary libraries
import streamlit as st
from src.chemprop1_helpers import *
import pandas as pd
import numpy as np
import plotly.express as px
import re

st.markdown("### Under construction")

# Streamlit Interface
st.markdown('## ChemProp2 Scoring Tool')

if 'nw' in st.session_state and 'chemprop_ft' in st.session_state and 'chemprop_md' in st.session_state:

    network_df = st.session_state['nw'].copy()
    features_df = st.session_state['chemprop_ft'].copy()
    metadata_df = st.session_state['chemprop_md'].copy()

    # Check if DataFrames are non-empty before processing
    if not network_df.empty and not features_df.empty and not metadata_df.empty:
        
        node_names = list(features_df.index)
    
        st.write("Running ChemProp2 Scoring...")
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


st.markdown('## ChemProp2 Scores Plot')

if 'ChemProp2_scores' in st.session_state and not st.session_state['ChemProp2_scores'].empty:
    scores_df = st.session_state['ChemProp2_scores'].copy()
    st.dataframe(scores_df)
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
