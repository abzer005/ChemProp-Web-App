# Import necessary libraries
import streamlit as st
from src.chemprop1_helpers import *
import pandas as pd

st.markdown("### ChemProp1 Analysis")

st.markdown("""
**Instructions:** 
For ChemProp1 analysis, select an attribute from your dataset that contains exactly two distinct sequential points (e.g., two different time points or two spatial points). 
This analysis is based on the assumption that the abundance of an initial compound decreases over time or across space, while the abundance of a new compound increases.
""")

# Check if metadata is loaded
if 'md_for_analysis' in st.session_state and not st.session_state['md_for_analysis'].empty:
    
    chemprop_md = st.session_state['md_for_analysis']

    df = inside_levels(chemprop_md)
    mask = df.apply(lambda row: len(row['LEVELS']) == 0, axis=1)
    df = df[~mask]
    
    # Display the overview dataframe and checkboxes
    st.markdown("### Select Metadata Columns for ChemProp1 Calculation")

    st.markdown("""
                - ChemProp1: Select the metadata column that contains the timepoint information. This column must be present in the metadata for ChemProp1 calculation and should contains ONLY 2 TIMEPOINTS.
                - Treatments: Select the metadata column that contains the treatment conditions (e.g., control and treatment groups). This category is optional.
                - Replicates: Select the metadata column that contains the replicate information. This category is optional and will be used for plotting intensity trends.
                """)
    
    st.markdown("""
                <div style="color: #05577C; font-size: 20px;">
                Select <b>only one metadata attribute</b> per checkbox column.
                </div>""", 
                unsafe_allow_html=True)

    # Add a new column for checkboxes (initialized to False)
    df['ChemProp1'] = False
    df['Treatments'] = False
    df['Replicates'] = False

    # Display an editable DataFrame using st.data_editor
    edited_df = st.data_editor(df, num_rows="fixed", 
                               use_container_width=True, 
                               hide_index = True,
                               disabled=("ATTRIBUTES", "LEVELS", "COUNTS"),
                               column_order=("ChemProp1", "Treatments" , "Replicates", "ATTRIBUTES", "LEVELS", "COUNTS")
                               )

    # Collect the rows where checkboxes are selected (True)
    chemprop_row = edited_df[edited_df['ChemProp1'] == True]['ATTRIBUTES'].tolist()
    treatment_row = edited_df[edited_df['Treatments'] == True]['ATTRIBUTES'].tolist()
    rep_row = edited_df[edited_df['Replicates'] == True]['ATTRIBUTES'].tolist()

    # Show the selected attributes for Treatments
    if chemprop_row:
         # Extract unique values from the selected attributes in `chemprop_md` dataframe
        unique_values = chemprop_md[chemprop_row].nunique().tolist()  # Number of unique values in each selected column

        # Display appropriate message based on unique values count
        if all(count < 3 for count in unique_values):  # Checking if all selected columns have fewer than 3 unique values
            st.success(
                f"You have only {', '.join(map(str, unique_values))} unique timepoints across selected attributes. "
                f"Hence, using ChemProp1 is appropriate for this data.",
                icon="✅"
            )
        else:
            st.error(
                "You cannot run ChemProp1 for attributes with more than 2 unique timepoints. "
                "Consider running ChemProp2, or filter the data further to have only 2 timepoints.",
                icon="🚨"
            )
        st.success(f"The column selected for calculating ChemProp1: {', '.join(chemprop_row)}")
    else:
        st.warning("No column selected for calculating ChemProp1 scores.")


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
    time_md = chemprop_md[chemprop_row].applymap(strip_non_numeric).apply(pd.to_numeric, errors='coerce')

    # Conditional creation of subset_chemprop_md
    if treatment_row and rep_row:
        subset_chemprop_md = pd.concat([time_md, treatment_md, rep_md], axis=1, join='outer')
    elif treatment_row:
        subset_chemprop_md = pd.concat([time_md, treatment_md], axis=1, join='outer')
    elif rep_row:
        subset_chemprop_md = pd.concat([time_md, rep_md], axis=1, join='outer')
    else:
        subset_chemprop_md = time_md
        
    st.markdown("### ChemProp1 Metadata with Selected Columns")
    st.dataframe(subset_chemprop_md)

    final_ft = st.session_state['ft_for_analysis'].copy()
    st.dataframe(final_ft)
    
    st.session_state['chemprop1_ft'] = final_ft
    st.session_state['chemprop1_subset_md'] = subset_chemprop_md
    st.session_state['chemprop1_md'] = time_md
else:
    st.warning("Metadata not loaded in the session state.")

################################## FILTER THE DATA ##########################################################
st.markdown("### Filter the Data")
subset_md = st.checkbox("Do you want to subset the data further?", False)

if subset_md:
    # Markdown for initial instruction
    st.markdown("Select samples based on their conditions.")

    # Container for the first filter - specifically for `chemprop_row` column
    with st.container():
        st.write("#### Initial Timepoint Filter")
        c1, c2 = st.columns(2)

        # First filter condition defaults to the timepoint column (`chemprop_row`)
        # This will be disabled so users cannot change it
        attribute_column = c1.selectbox(
            "Select the column for categorization (Timepoints)", 
            options=chemprop_row, 
            disabled=True
        )

        # Multi-select for exactly two timepoints
        subset_group = c2.multiselect(
            "Select exactly TWO timepoints for filtering", 
            options=sorted(subset_chemprop_md[chemprop_row].dropna().astype(str).stack().unique()),
            max_selections=2
        )

        # Apply the first filter if two timepoints are selected
        if len(subset_group) == 2:

            # Filter indices based on the selected timepoints
            # Ensure `timepoint_column` is correctly set to the column name from `chemprop_row`
            timepoint_column = chemprop_row if isinstance(chemprop_row, str) else chemprop_row[0]
            subset_group = list(map(str, subset_group))  # Convert subset_group to strings if needed for matching
            subgroup_indices = subset_chemprop_md[subset_chemprop_md[timepoint_column].astype(str).isin(subset_group)].index

            # Filter feature table and metadata based on the selected indices
            subset_ft = st.session_state['ft_for_analysis'].loc[:, subgroup_indices]
            subset_md = subset_chemprop_md.loc[subgroup_indices]

            # Display the filtered data immediately
            with st.expander(f"Filtered Group {subset_ft.shape}"):
                st.dataframe(subset_ft)
                st.dataframe(subset_md)

            # Save to session state for later use
            st.session_state['chemprop1_ft'] = subset_ft
            st.session_state['chemprop1_subset_md'] = subset_md
            st.session_state['chemprop1_md'] = subset_md[chemprop_row]

            # Option to add another filter
            add_filter = st.checkbox("Would you like to add another filter?", False)
            if add_filter:
                st.markdown("#### Additional Filter")
                with st.container():
                    c1, c2 = st.columns(2)

                    # Allow the user to select any column for further filtering
                    additional_column = c1.selectbox(
                        "Select the column for categorization",
                        options=st.session_state['chemprop1_subset_md'].columns,
                        key="additional_column"
                    )

                    # Multi-select for categories in the selected column
                    additional_group = c2.multiselect(
                        "Select categories for additional filtering",
                        options=sorted(st.session_state['chemprop1_subset_md'][additional_column].dropna().unique()),
                        key="additional_group"
                    )

                    # Apply the additional filter if categories are selected
                    if additional_group:
                       # timepoint_column = chemprop_row if isinstance(chemprop_row, str) else chemprop_row[0]
                                               
                       # subgroup_indices = subset_chemprop_md[subset_chemprop_md[timepoint_column].astype(str).isin(subset_group)].index

                        additional_group = list(map(str, additional_group))  # Convert the group to strings if needed for matching
                        additional_indices = st.session_state['chemprop1_subset_md'][st.session_state['chemprop1_subset_md'][additional_column].astype(str).isin(additional_group)].index
                        
                        # Update the feature table and metadata based on additional filtering
                        final_ft = subset_ft.loc[:, additional_indices]
                        final_md = subset_md.loc[additional_indices]

                        # Display the final filtered data
                        with st.expander(f"Final Filtered Group {final_ft.shape}"):
                            st.dataframe(final_ft)
                            st.dataframe(final_md)

                        # Update session state with the final filtered tables
                        st.session_state['chemprop1_ft'] = final_ft
                        st.session_state['chemprop1_subset_md'] = final_md
                        st.session_state['chemprop1_md'] = final_md[chemprop_row]

# Streamlit Interface
st.markdown('## ChemProp1 Scores')

if 'nw' in st.session_state and 'chemprop1_ft' in st.session_state and 'chemprop1_md' in st.session_state:

    network_df = st.session_state['nw'].copy()
    features_df = st.session_state['chemprop1_ft'].copy()
    metadata_df = st.session_state['chemprop1_md'].copy()

    # Check if DataFrames are non-empty before processing
    if not network_df.empty and not features_df.empty and not metadata_df.empty:

        # Add a button to trigger ChemProp2 Scoring
        if "run_chemprop1" not in st.session_state:
            st.session_state.run_chemprop1 = False

        if st.button("Run ChemProp1") or st.session_state.run_chemprop1:
            st.session_state.run_chemprop1 = True

            # Fix the function and the way result_df will look like
            chemprop1_result_df = run_chemprop1_analysis(network_df, features_df, metadata_df)

            # Get the number of columns in network_df
            network_cols = network_df.shape[1]

            # Define the range of columns in result_df you want to check (from 6th column to 11th column)
            extra_cols_range = chemprop1_result_df.columns[network_cols:]  # From ncol(network_df) + 1 to ncol(result_df)

            # Drop rows where all values in the selected columns are NaN
            chemprop1_result_df = chemprop1_result_df.dropna(subset=extra_cols_range, how='all')

            # Check if 'an_gnps' is present in the session state and is not empty
            if 'an_gnps' in st.session_state and not st.session_state['an_gnps'].empty:
    
                # Get the DataFrames from session state
                gnps_df = st.session_state['an_gnps'].copy()

                # Call the function to add the 'ID1_name' and 'ID2_name' columns
                updated_chemprop1_df = add_names_to_chemprop(chemprop1_result_df, gnps_df)

                # Save the updated DataFrame back to session state
                st.session_state['ChemProp1_scores'] = updated_chemprop1_df

               # Update result_df to reflect the changes
                chemprop1_result_df = updated_chemprop1_df

            # Display the (potentially updated) result_df
            st.write("ChemProp1 Scoring Results:")
            st.dataframe(chemprop1_result_df)
            st.session_state['ChemProp1_scores'] = chemprop1_result_df

            # Display only the row counts for both DataFrames
            st.markdown(
                f"""
                <p style="color:red;">Number of edges in the original edge table: <b>{network_df.shape[0]}</b></p>
                <p style="color:red;">Number of edges in the ChemProp1 edge table: <b>{chemprop1_result_df.shape[0]}</b></p>
                """,
                unsafe_allow_html=True)
            
            if network_df.shape[0] != chemprop1_result_df.shape[0]:
                st.warning('The reduced number of edges might be due to the removal of blank entries.')

            # Allow users to download the results
            # Let the user input the filename (default value provided)
            user_filename = st.text_input("Enter the filename for the CSV and press Enter to apply the name:", 
                                          value="chemprop1_scores_results.csv")

            # Ensure the file ends with '.csv'
            if not user_filename.endswith('.csv'):
                user_filename += '.csv'
                
            # Convert result_df to CSV and encode it
            chemprop1_result = chemprop1_result_df.to_csv(index=False).encode('utf-8')

            # Create three columns in the UI for buttons
            col1, col2, col3, col4 = st.columns([1.1, 1, 1, 4.9]) # to make the buttons lie close to each other

            # Column 1: Download CSV
            with col1:
                # Provide a download button for CSV
                st.download_button(
                    label="Download Results as CSV",
                    data=chemprop1_result,  
                    file_name=user_filename,
                    mime='text/csv'
                )

            # Column 2: Download GraphML
            with col2:
                # Button to trigger GraphML generation
                if __name__ == '__main__':
                    generate_graphml_zip_chemprop1()  # Function to generate and download GraphML file

            # Column 3: Visualize Network
            with col3:
                # Button to visualize the network
                if st.button("Visualize Network in GNPS2"):
                    if 'ChemProp1_scores' in st.session_state:
                        st.write('This functionality will be added soon')
                         # Call function to visualize the network


            # st.markdown('### Scores Plot Visualization')

            #Getting the scores plot:
            scores1_df = st.session_state['ChemProp1_scores'].copy()
            scores1_df['CLUSTER IDs'] = scores1_df['CLUSTERID1'].astype(str) + '-' + scores1_df['CLUSTERID2'].astype(str)

            # Plot a scatter plot using Plotly
            fig1 = px.scatter(
                scores1_df, 
                x='DeltaMZ', 
                y='ChemProp1', 
                hover_name='CLUSTER IDs',  # Show the 'Labels' column when hovering over a point
                title='Scatter Plot: DeltaMZ vs ChemProp1 scores',
                labels={'DeltaMZ': 'Delta M/Z', 'ChemProp1': 'ChemProp1 scores'}
                )

            # Show the plot in Streamlit
            st.plotly_chart(fig1)

            # Check if the plot exists (based on session state or similar condition)
            if 'ChemProp1_scores' in st.session_state and st.session_state.run_chemprop1:
                show_filters = st.checkbox('Show Filters')

                # If the checkbox is checked, show additional options or success message
                if show_filters:

                    filter_option_chemprop1 = st.radio(
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
                    if filter_option_chemprop1 == 'Filter by Score':
                        with c1:
                            score1_min = st.number_input(f"Min Score:", value=float(scores1_df['ChemProp1'].min()))
                        with c2:
                            score1_max = st.number_input(f"Max Score:", value=float(scores1_df['ChemProp1'].max()))
                        
                        select_df_chemprop1 = scores1_df[(scores1_df['ChemProp1'] >= score1_min) & (scores1_df['ChemProp1'] <= score1_max)]
                        
                        with c3:
                            # Create concatenated column for dropdown display with chemprop1 rounded to 3 decimal places
                            select_df_chemprop1['Dropdown_Display'] = "ID: " + select_df_chemprop1['CLUSTER IDs'].astype(str) + ", Score: " + select_df_chemprop1['ChemProp1'].round(3).astype(str)

                            # Create a dropdown with the concatenated names
                            selected_option_chemprop1_chemprop1 = st.selectbox("Select the edge to view plots", select_df_chemprop1['Dropdown_Display'].tolist())

                            if selected_option_chemprop1_chemprop1 in select_df_chemprop1['Dropdown_Display'].values:
                                selected_row_data = select_df_chemprop1.loc[select_df_chemprop1['Dropdown_Display'] == selected_option_chemprop1_chemprop1]
                                selected_row_data = selected_row_data.drop(columns=['Dropdown_Display'])
                                filtered_df = pd.DataFrame(selected_row_data)

                    elif filter_option_chemprop1 == 'Filter by M/Z Range':
                        with c1:
                            mz_min = st.number_input("Min DeltaMZ:", value=float(scores1_df['DeltaMZ'].min()))
                        with c2:
                            mz_max = st.number_input("Max DeltaMZ:", value=float(scores1_df['DeltaMZ'].max()))

                        select_df_chemprop1 = scores1_df[
                             (scores1_df['DeltaMZ'] >= mz_min) & (scores1_df['DeltaMZ'] <= mz_max)
                         ]
                        
                        with c3:
                            # Create concatenated column for dropdown display
                            select_df_chemprop1['Dropdown_Display'] = "ID: " + select_df_chemprop1['CLUSTER IDs'].astype(str) + ", DeltaMZ: " + select_df_chemprop1['DeltaMZ'].round(3).astype(str)

                            # Create a dropdown with the concatenated names
                            selected_option_chemprop1 = st.selectbox("Select the edge to view plots", select_df_chemprop1['Dropdown_Display'].tolist())

                            if selected_option_chemprop1 in select_df_chemprop1['Dropdown_Display'].values:
                                selected_row_data = select_df_chemprop1.loc[select_df_chemprop1['Dropdown_Display'] == selected_option_chemprop1]
                                selected_row_data = selected_row_data.drop(columns=['Dropdown_Display'])
                                filtered_df = pd.DataFrame(selected_row_data)
                        
                    elif filter_option_chemprop1 == 'Filter by Name':
                        # Check if ID1_Name and ID2_Name columns exist
                        if 'ID1_name' in scores1_df.columns and 'ID2_name' in scores1_df.columns:
                            with c1:
                                cluster_name = st.text_input("Enter text:", "")
                                select_df_chemprop1 = scores1_df.copy()
                            
                                if cluster_name:
                                        # Filter rows where the cluster_name is found in either ID1_Name or ID2_Name
                                        select_df_chemprop1 = scores1_df[
                                            scores1_df['ID1_name'].str.contains(cluster_name, case=False, na=False) |
                                            scores1_df['ID2_name'].str.contains(cluster_name, case=False, na=False)
                                        ]

                            with c2:
                                # Create concatenated column for dropdown display
                                select_df_chemprop1['Dropdown_Display'] = "ID: " + select_df_chemprop1['CLUSTER IDs'].astype(str)

                                # Create a dropdown with the concatenated names
                                selected_option_chemprop1 = st.selectbox("Select the edge to view plots", select_df_chemprop1['Dropdown_Display'].tolist())

                                if selected_option_chemprop1 in select_df_chemprop1['Dropdown_Display'].values:
                                    selected_row_data = select_df_chemprop1.loc[select_df_chemprop1['Dropdown_Display'] == selected_option_chemprop1]
                                    selected_row_data = selected_row_data.drop(columns=['Dropdown_Display'])
                                    filtered_df = pd.DataFrame(selected_row_data)     

                        else:
                            # Show a message if the names columns are not available
                            st.warning("You don't have names for the cluster IDs in the table.")

                    elif filter_option_chemprop1 == 'Filter by Cluster ID':     
                        
                        # Dropdown list for selecting a combined CLUSTER_ID
                        selected_cluster_id = st.selectbox("Select the edge to view plots:", options=scores1_df['CLUSTER IDs'].unique())

                        # Filter the DataFrame based on the selected Cluster ID
                        filtered_df = scores1_df[scores1_df['CLUSTER IDs'] == selected_cluster_id]

                        # Extract the actual CLUSTERID1 and CLUSTERID2 values from the selected_cluster_id
                        selected_cluster_values = selected_cluster_id.split('-')
                        cluster_id1 = int(selected_cluster_values[0])
                        cluster_id2 = int(selected_cluster_values[1])

                    # Replot the filtered data
                    st.dataframe(filtered_df)
                    
                    # Create two rows with columns, setting appropriate proportions
                    with st.container():
                        c1, c2 = st.columns(2)  
                        
                        with c1:
                            if not filtered_df.empty:
                                selected_row = filtered_df.iloc[0]  # Select the first (or only) row
                                fig = plot_intensity_trends_single_row_chemprop1(selected_row, features_df, st.session_state['chemprop1_md'])
                
                                st.plotly_chart(fig, use_container_width=True)

                        with c2:
                            try:
                                # Check if ChemProp1_scores is available in session state
                                if 'ChemProp1_scores' in st.session_state and not filtered_df.empty:
                                    df = st.session_state['ChemProp1_scores'].copy()
                                    comp_index = filtered_df['ComponentIndex'].iloc[0]
                                    plot_df = df[df['ComponentIndex'] == comp_index]

                                    # List of allowed columns for edge labels
                                    allowed_columns = [
                                        'ComponentIndex', 'Cosine', 'DeltaMZ', 
                                        'ChemProp1', 'Sign_ChemProp1', 'abs_ChemProp1'
                                        ]
                                
                                    available_columns = [col for col in plot_df.columns if col in allowed_columns]

                                    # Let the user select which column to use as edge labels
                                    edge_label_column = st.selectbox("Select column for edge labels:", 
                                                                     options=available_columns,
                                                                     help="To save the network image as PNG, right-click on the empty space and select 'Save image as'."
                                    )
                                    # Generate the graph using the selected column, id1, and id2
                                    nodes, edges = generate_graph_from_df_chemprop1(plot_df, filtered_df, edge_label_column)

                                    # Define graph configuration
                                    config = Config(width=800, 
                                                    height=600, 
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
                                    st.write("No ChemProp1 scores data available.")

                            except ModuleNotFoundError:
                                st.error("This page requires the `pygraphviz` package, which is not available in the Windows app.")