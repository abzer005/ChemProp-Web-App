# Import necessary libraries
import streamlit as st
from src.chemprop1_helpers import *
import pandas as pd
import streamlit_agraph

page_setup()

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
    df['ChemProp1 (Timepoint)'] = False
    df['Treatments'] = False
    df['Replicates'] = False

    # Display an editable DataFrame using st.data_editor
    edited_df = st.data_editor(df, num_rows="fixed", 
                               use_container_width=True, 
                               hide_index = True,
                               disabled=("ATTRIBUTES", "LEVELS", "COUNTS"),
                               column_order=("ChemProp1 (Timepoint)", "Treatments" , "Replicates", "ATTRIBUTES", "LEVELS", "COUNTS")
                               )

    # Collect the rows where checkboxes are selected (True)
    chemprop_row = edited_df[edited_df['ChemProp1 (Timepoint)'] == True]['ATTRIBUTES'].tolist()
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

    final_ft = st.session_state['ft_for_analysis'].copy()
    st.session_state['chemprop1_ft'] = final_ft
    st.session_state['chemprop1_subset_md'] = subset_chemprop_md
    st.session_state['chemprop1_md'] = time_md

    tab_info= [
        ("chemprop1_subset_md", "Filtered Metadata"),
        ("chemprop1_ft", "Filtered Feature Table"),
    ]

    display_tables_in_tabs(tab_info)

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
            
            # Save to session state for later use
            st.session_state['chemprop1_ft'] = subset_ft
            st.session_state['chemprop1_subset_md'] = subset_md
            st.session_state['chemprop1_md'] = subset_md[chemprop_row]

            filtered_tab_info= [
                ("chemprop1_subset_md", "Filtered Metadata"),
                ("chemprop1_ft", "Filtered Feature Table"),
            ]

            # Display the filtered data immediately
            with st.expander(f"Filtered Tables"):
                display_tables_in_tabs(filtered_tab_info)

            # Option to add another filter
            add_filter = st.checkbox("Would you like to add another filter?", False)
            if add_filter:
                with st.container():
                    c1, c2 = st.columns(2)

                    # Allow the user to select any column for further filtering
                    additional_column = c1.selectbox(
                        "Select the column for categorization (Ex: Treatments, Replicates)",
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
                
                        # Update session state with the final filtered tables
                        st.session_state['chemprop1_ft'] = final_ft
                        st.session_state['chemprop1_subset_md'] = final_md
                        st.session_state['chemprop1_md'] = final_md[chemprop_row]

                        # Display the final filtered data
                        final_filtered_tabs= [
                            ("chemprop1_subset_md", "Filtered Metadata"),
                            ("chemprop1_ft", "Filtered Feature Table"),
                        ]

                        # Display the filtered data immediately
                        with st.expander(f"Final Filtered Tables"):
                            display_tables_in_tabs(final_filtered_tabs)

# Streamlit Interface
st.markdown('## Calculate ChemProp1 Scores')

dfs = get_inputs_from_state()

if dfs is None:
    st.info("Upload Edge Table + Feature Table + Metadata to run ChemProp1.")
    st.stop()

network_df, features_df, metadata_df = dfs

run_flag, mode = chemprop1_controls()
if not run_flag:
    st.stop()

if "gnps_task_id" in st.session_state and st.session_state["gnps_task_id"]:
    st.session_state["gnps_task_id"] = str(st.session_state["gnps_task_id"]).strip()

if "last_chemprop1_mode" not in st.session_state:
    st.session_state["last_chemprop1_mode"] = mode

if st.session_state["last_chemprop1_mode"] != mode:
    st.session_state["show_filters_chemprop1"] = False  # reset checkbox
    st.session_state["last_chemprop1_mode"] = mode
    st.rerun()

if mode == "Provided Edge Table":
        chemprop1_result_df = run_chemprop1_pipeline(network_df, features_df, metadata_df)
       
        if chemprop1_result_df is None or chemprop1_result_df.empty:
            st.warning("ChemProp1 returned no results after filtering.")
            st.stop()

        render_results_summary(network_df, chemprop1_result_df)
        render_download_buttons(chemprop1_result_df)
        render_scores_plot(chemprop1_result_df)

        render_filters_and_plots(
            scores_df=chemprop1_result_df,
            features_df=features_df,
            metadata_df=metadata_df,
        )
        

elif mode == "Cascade edges":

    # Extract feature IDs from feature table
    feature_ids = (
        features_df.index
        .dropna()
        .astype(int)
        .sort_values()
        .unique()
        .tolist()
    )

    source_id = st.selectbox(
        "Select Feature ID (source)",
        options=feature_ids,
        index=0  
    )

    cascade_df = build_cascade_edges(
        edges=network_df,
        source_id=int(source_id),
        max_neighbor=15,
        col1="CLUSTERID1",
        col2="CLUSTERID2",
        component_col="ComponentIndex",
    )

    if cascade_df.empty:
        st.warning(f"No cascade edges found for Feature ID {source_id}.")
        st.stop()

    with st.expander("Cascade edge table (source → neighbors)"):
        st.caption(
            f"{cascade_df.shape[0]} rows × {cascade_df.shape[1]} columns · "
            "Edges ordered by neighbor level"
        )
        st.dataframe(cascade_df, hide_index=True, use_container_width=True)

    chemprop1_result_df = run_chemprop1_pipeline(cascade_df, features_df, metadata_df)
    
    render_results_summary(cascade_df, chemprop1_result_df)
    render_download_buttons(chemprop1_result_df)
    render_scores_plot(chemprop1_result_df)

    render_filters_and_plots(
        scores_df=chemprop1_result_df,
        features_df=features_df,
        metadata_df=metadata_df,
    )

elif mode == "User defined edge":
    
    colA, colB = st.columns(2)

    # Use feature IDs from the feature table index as choices (safer than free typing)
    ft_ids = features_df.index.astype(str).tolist()

    with colA:
        id1_str = st.selectbox("CLUSTERID1 (Source)", options=ft_ids, index=0)
    with colB:
        id2_str = st.selectbox("CLUSTERID2 (Target)", options=ft_ids, index=min(1, len(ft_ids)-1))

    # Convert back to same dtype as features_df index if needed
    try:
        id1 = int(id1_str)
        id2 = int(id2_str)
    except Exception:
        id1, id2 = id1_str, id2_str

    if id1 == id2:
        st.warning("Please choose two different feature IDs.")
        st.stop()
    
    user_edge_df = pd.DataFrame(
        [{
            "CLUSTERID1": id1,
            "CLUSTERID2": id2,
            "ComponentIndex": 0,  # dummy
            "Cosine": np.nan,      # optional placeholder
        }]
    )

    chemprop1_result_df = run_chemprop1_pipeline(user_edge_df, features_df, metadata_df)

    if chemprop1_result_df is None or chemprop1_result_df.empty:
        st.warning("ChemProp1 returned no results for the selected edge (feature missing after filtering?).")
        st.stop()

    render_results_summary(user_edge_df, chemprop1_result_df)
   
    render_edge_detail_plots(
        filtered_df=chemprop1_result_df,
        all_scores_df=chemprop1_result_df,
        features_df=features_df,
        metadata_df=metadata_df,
    )
    if (st.session_state.get("gnps_task_id") or "").strip():
        render_spectra_modifinder(chemprop1_result_df, st.session_state.get("an_gnps"))

