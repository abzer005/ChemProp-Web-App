# Import necessary libraries
import streamlit as st
from src.chemprop1_helpers import *
from src.chemprop2_helpers import *
import pandas as pd
import re
import plotly.express as px


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
                - Replicates: Select the metadata column that contains the replicate information. This category is optional and will be used for plotting intensity trends.
                """)
    
    st.markdown("""
                <div style="color: #05577C; font-size: 20px;">
                Select <b>only one metadata attribute</b> per checkbox column.
                </div>""", 
                unsafe_allow_html=True)

    # Add a new column for checkboxes (initialized to False)
    df['ChemProp2'] = False
    df['Treatments'] = False
    df['Replicates'] = False

    # Display an editable DataFrame using st.data_editor
    edited_df = st.data_editor(df, num_rows="fixed", 
                               use_container_width=True, 
                               hide_index = True,
                               disabled=("ATTRIBUTES", "LEVELS", "COUNTS"),
                               column_order=("ChemProp2", "Treatments" , "Replicates", "ATTRIBUTES", "LEVELS", "COUNTS")
                               )

    # Collect the rows where checkboxes are selected (True)
    chemprop2_row = edited_df[edited_df['ChemProp2'] == True]['ATTRIBUTES'].tolist()
    treatment_row = edited_df[edited_df['Treatments'] == True]['ATTRIBUTES'].tolist()
    rep_row = edited_df[edited_df['Replicates'] == True]['ATTRIBUTES'].tolist()

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

    # Show the selected attributes for Replicates
    if rep_row:
        st.success(f"The column selected for Replicate conditions: {', '.join(rep_row)}")
    else:
        st.warning("No column selected for Replicates.")
    
    # Subset the DataFrames based on the user's selection
    treatment_md = chemprop_md[treatment_row] if treatment_row else pd.DataFrame()
    rep_md = chemprop_md[rep_row] if rep_row else pd.DataFrame()
    time_md = chemprop_md[chemprop2_row]

    # Apply the function to the subset DataFrame
    time_md = time_md.applymap(strip_non_numeric)

    # Convert the stripped values to numeric types (float)
    time_md = time_md.apply(pd.to_numeric, errors='coerce')

    # Conditional creation of subset_chemprop_md
    if treatment_row and rep_row:
        subset_chemprop_md = pd.concat([time_md, treatment_md, rep_md], axis=1, join='outer')
    elif treatment_row:
        subset_chemprop_md = pd.concat([time_md, treatment_md], axis=1, join='outer')
    elif rep_row:
        subset_chemprop_md = pd.concat([time_md, rep_md], axis=1, join='outer')
    else:
        subset_chemprop_md = time_md
        
    final_ft = st.session_state['ft_for_analysis'].copy()
    
    with st.expander(f"ChemProp2 Data with Selected Columns"):
                            st.dataframe(subset_chemprop_md)
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
        if "run_chemprop2" not in st.session_state:
            st.session_state.run_chemprop2 = False

        if st.button("Run ChemProp2") or st.session_state.run_chemprop2:
            st.session_state.run_chemprop2 = True
            node_names = list(features_df.index)
    
            result_df = ChemProp2_scoring(network_df, features_df, node_names, metadata_df)
        
            st.write("ChemProp2 Scoring Results:")

            # Get the number of columns in network_df
            network_cols = network_df.shape[1]

            # Define the range of columns in result_df you want to check (from 6th column to 11th column)
            extra_cols_range = result_df.columns[network_cols:]  # From ncol(network_df) + 1 to ncol(result_df)

            # Drop rows where all values in the selected columns are NaN
            result_df = result_df.dropna(subset=extra_cols_range, how='all')

            # Check if 'an_gnps' is present in the session state and is not empty
            if 'an_gnps' in st.session_state and st.session_state['an_gnps'] is not None and not st.session_state['an_gnps'].empty:

                # Get the DataFrames from session state
                gnps_df = st.session_state['an_gnps'].copy()

                # Call the function to add the 'ID1_name' and 'ID2_name' columns
                updated_chemprop_df = add_names_to_chemprop(result_df, gnps_df)

                # Save the updated DataFrame back to session state
                st.session_state['ChemProp2_scores'] = updated_chemprop_df

               # Update result_df to reflect the changes
                result_df = updated_chemprop_df

                st.dataframe(result_df)
                st.session_state['ChemProp2_scores'] = result_df
            else:
                # Display the (potentially updated) result_df
                st.dataframe(result_df)
                st.session_state['ChemProp2_scores'] = result_df

            # Display only the row counts for both DataFrames
            st.markdown(
                f"""
                <p style="color:red;">Number of edges in the original edge table: <b>{network_df.shape[0]}</b></p>
                <p style="color:red;">Number of edges in the ChemProp2 edge table: <b>{result_df.shape[0]}</b></p>
                """,
                unsafe_allow_html=True)
            
            if network_df.shape[0] != result_df.shape[0]:
                st.warning('The reduced number of edges might be due to the removal of blank entries.')
          
            user_filename = st.text_input("**Enter the filename for the CSV and press Enter to apply the name:**", 
                                          value="chemprop2_scores_results.csv")
                    
            # Check if the user has entered a filename
            if not user_filename.endswith('.csv'):
                user_filename += '.csv'

            # Convert result_df to CSV and encode it
            chemprop2_result = result_df.to_csv(index=False).encode('utf-8')

            # Create three columns in the UI for buttons
            col1, col2, col3, col4 = st.columns([1.1, 1.1, 1.1, 4.5]) # to make the buttons lie close to each other

            with col1:
                st.download_button(
                    label="Download result as CSV",
                    data=chemprop2_result,
                    file_name=user_filename,
                    mime='text/csv'
                )

            # Column 2: Download GraphML
            with col2:
                # Button to trigger GraphML generation
                if __name__ == '__main__':
                    generate_graphml_zip_chemprop2()

            # Column 3: Visualize Network
            with col3:
                # Button to visualize the network
                if st.button("Visualize Network in GNPS2"):
                    if 'ChemProp2_scores' in st.session_state:
                        st.write('This functionality will be added soon')
                         # Call function to visualize the network
            
            ##########################

            st.markdown('### False Discovery Rate')

            if st.button("Apply FDR"):
                if 'ChemProp2_scores' in st.session_state:
                    st.write('Select the positive and negative cutoffs for the ChemProp2 scores based on your FDR-curve')
                    overall_fdr_table, fig_histogram, fig_fdr = calculate_fdr(network_df, features_df, metadata_df, score_range=(-1, 1), bin_size=0.001)

                     # Sample and blank selection interface
                    c1, c2 = st.columns(2)
                    c1.plotly_chart(fig_histogram)
                    c2.plotly_chart(fig_fdr)
                    
            ##########################

            st.markdown('### Scores Plot Visualization')

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

            # Check if the plot exists (based on session state or similar condition)
            if 'ChemProp2_scores' in st.session_state and st.session_state.run_chemprop2:
                show_filters = st.checkbox('Show Filters')

                # If the checkbox is checked, show additional options or success message
                if show_filters:

                    filter_option = st.radio(
                            "Filter type:",
                            ('Filter by Score', 
                             'Filter by M/Z Range',
                             'Filter by Name',
                             'Filter by Cluster ID'),
                             horizontal = True
                        )

                    # Create three columns
                    c1, c2, c3 = st.columns(3)

                    # Show the appropriate input fields based on the selected filter option
                    if filter_option == 'Filter by Score':
                        with c1:
                            score_min = st.number_input(f"Min {y_axis_option}:", value=float(scores_df['ChemProp2'].min()))
                        with c2:
                            score_max = st.number_input(f"Max {y_axis_option}:", value=float(scores_df['ChemProp2'].max()))
                        
                        select_df = scores_df[(scores_df[y_axis_option] >= score_min) & (scores_df[y_axis_option] <= score_max)]
                        
                        with c3:
                            # Create concatenated column for dropdown display with y_axis_option rounded to 3 decimal places
                            select_df['Dropdown_Display'] = "ID: " + select_df['CLUSTER IDs'].astype(str) + ", Score: " + select_df[y_axis_option].round(3).astype(str)

                            # Create a dropdown with the concatenated names
                            selected_option = st.selectbox("Select the edge to view plots", select_df['Dropdown_Display'].tolist())

                            if selected_option in select_df['Dropdown_Display'].values:
                                selected_row_data = select_df.loc[select_df['Dropdown_Display'] == selected_option]
                                selected_row_data = selected_row_data.drop(columns=['Dropdown_Display'])
                                filtered_df = pd.DataFrame(selected_row_data)

                    elif filter_option == 'Filter by M/Z Range':
                        with c1:
                            mz_min = st.number_input("Min DeltaMZ:", value=float(scores_df['DeltaMZ'].min()))
                        with c2:
                            mz_max = st.number_input("Max DeltaMZ:", value=float(scores_df['DeltaMZ'].max()))

                        select_df = scores_df[
                             (scores_df['DeltaMZ'] >= mz_min) & (scores_df['DeltaMZ'] <= mz_max)
                         ]
                        
                        with c3:
                            # Create concatenated column for dropdown display
                            select_df['Dropdown_Display'] = "ID: " + select_df['CLUSTER IDs'].astype(str) + ", DeltaMZ: " + select_df['DeltaMZ'].round(32).astype(str)

                            # Create a dropdown with the concatenated names
                            selected_option = st.selectbox("Select the edge to view plots", select_df['Dropdown_Display'].tolist())

                            if selected_option in select_df['Dropdown_Display'].values:
                                selected_row_data = select_df.loc[select_df['Dropdown_Display'] == selected_option]
                                selected_row_data = selected_row_data.drop(columns=['Dropdown_Display'])
                                filtered_df = pd.DataFrame(selected_row_data)
                        
                    elif filter_option == 'Filter by Name':
                        # Check if ID1_Name and ID2_Name columns exist
                        if 'ID1_name' in scores_df.columns and 'ID2_name' in scores_df.columns:
                            with c1:
                                cluster_name = st.text_input("Enter text:", "")
                                select_df = scores_df.copy()
                            
                                if cluster_name:
                                        # Filter rows where the cluster_name is found in either ID1_Name or ID2_Name
                                        select_df = scores_df[
                                            scores_df['ID1_name'].str.contains(cluster_name, case=False, na=False) |
                                            scores_df['ID2_name'].str.contains(cluster_name, case=False, na=False)
                                        ]

                            with c2:
                                # Create concatenated column for dropdown display
                                select_df['Dropdown_Display'] = "ID: " + select_df['CLUSTER IDs'].astype(str)

                                # Create a dropdown with the concatenated names
                                selected_option = st.selectbox("Select the edge to view plots", select_df['Dropdown_Display'].tolist())

                                if selected_option in select_df['Dropdown_Display'].values:
                                    selected_row_data = select_df.loc[select_df['Dropdown_Display'] == selected_option]
                                    selected_row_data = selected_row_data.drop(columns=['Dropdown_Display'])
                                    filtered_df = pd.DataFrame(selected_row_data)     

                        else:
                            # Show a message if the names columns are not available
                            st.warning("You don't have names for the cluster IDs in the table.")

                    elif filter_option == 'Filter by Cluster ID':     
                        
                        # Dropdown list for selecting a combined CLUSTER_ID
                        selected_cluster_id = st.selectbox("Select the edge to view plots:", options=scores_df['CLUSTER IDs'].unique())

                        # Filter the DataFrame based on the selected Cluster ID
                        filtered_df = scores_df[scores_df['CLUSTER IDs'] == selected_cluster_id]

                        # Extract the actual CLUSTERID1 and CLUSTERID2 values from the selected_cluster_id
                        selected_cluster_values = selected_cluster_id.split('-')
                        cluster_id1 = int(selected_cluster_values[0])
                        cluster_id2 = int(selected_cluster_values[1])

                    # Replot the filtered data
                    st.dataframe(filtered_df)
                    
                    # Create two rows with columns, setting appropriate proportions
                    with st.container():
                        c1, c2 = st.columns(2)  

                        # Row 1, Column 1: Plot intensity trend
                        with c1:
                            if not filtered_df.empty:
                                selected_row = filtered_df.iloc[0]  # Select the first (or only) row
                                fig = plot_intensity_trends_single_row_chemprop2(selected_row, features_df, st.session_state['chemprop_md'])
                                st.plotly_chart(fig, use_container_width=True)

                        # Row 1, Column 2
                        with c2:
                            try:
                                # Check if ChemProp2_scores is available in session state
                                if 'ChemProp2_scores' in st.session_state and not filtered_df.empty:
                                    df = st.session_state['ChemProp2_scores'].copy()
                                    comp_index = filtered_df['ComponentIndex'].iloc[0]
                                    plot_df = df[df['ComponentIndex'] == comp_index]

                                    # List of allowed columns for edge labels
                                    allowed_columns = [
                                        'ComponentIndex', 'Cosine', 'DeltaMZ', 'ChemProp2', 
                                        'ChemProp_spearman', 'ChemProp_log', 'ChemProp_sqrt', 
                                        'Sign_ChemProp2', 'abs_ChemProp2', 'abs_ChemProp_spearman', 
                                        'abs_ChemProp_log', 'abs_ChemProp_sqrt'
                                        ]
                                
                                    available_columns = [col for col in plot_df.columns if col in allowed_columns]
                                
                                    # Let the user select which column to use as edge labels
                                    edge_label_column = st.selectbox("Select column for edge labels:", 
                                                                     options=available_columns,                                                                   
                                                                     help="To save the network image as PNG, right-click on the empty space and select 'Save image as'."
                                    )

                                    # Generate the graph using the selected column, id1, and id2
                                    nodes, edges = generate_graph_from_df_chemprop2(plot_df, filtered_df, edge_label_column)

                                    # Define graph configuration
                                    config = Config(width=800, height=600, 
                                                    directed=True, 
                                                    nodeHighlightBehavior=True, 
                                                    highlightColor="#F7A7A6", 
                                                    collapsible=True, 
                                                    node={'labelProperty': 'label'}, 
                                                    link={'labelProperty': 'label', 
                                                          'renderLabel': True},
                                                    staticGraph=False)

                                    # Custom HTML legend with blue and red circles and text inside
                                    st.markdown("""
                                    <div style="display: flex; justify-content: center;">
                                        <div style="display: flex; align-items: center; margin-right: 20px;">
                                            <div style="width: 40px; height: 40px; border-radius: 50%; background-color: blue; color: white; display: flex; align-items: center; justify-content: center; font-size: 12px;">
                                                Source
                                            </div>
                                        </div>
                                        <div style="display: flex; align-items: center;">
                                            <div style="width: 40px; height: 40px; border-radius: 50%; background-color: red; color: white; display: flex; align-items: center; justify-content: center; font-size: 12px;">
                                                Target
                                            </div>
                                        </div>
                                    </div>
                                    """, unsafe_allow_html=True)

                                    # Display the graph in Streamlit
                                    agraph(nodes=nodes, edges=edges, config=config)
    
                                else:
                                    st.write("No ChemProp2 scores data available.")
                            
                            except ModuleNotFoundError:
                                st.error("This page requires the `pygraphviz` package, which is not available in the Windows app.")
                        




        