# Import necessary libraries
import streamlit as st
from src.chemprop1_helpers import *
from src.chemprop2_helpers import *
import pandas as pd
import re
import plotly.express as px

page_setup()

st.markdown("## ChemProp2 Analysis for Time Series Data")

with st.expander("ℹ️ ChemProp2 – How it works", expanded=False):
    st.markdown(
        """
        🔹 **Choose an attribute**  
        - Select a variable with **more than two timepoints**  
        - ChemProp2 Analyzes **pairs of features** and how they change **relative to each other**

        🔹 **Correlation logic**
        - 📈📈 **Move together** → positive correlation → **score ≈ 0**
        - 📈📉 **Move opposite** → anti-correlation → **directional score**

        🔹 **ChemProp2 Score (–1 to +1)**
        - 💪 **|Score|** → strength of potential biotransformation  
        - ➡️ **Positive (+)** → **Feature A → Feature B**  
        - ⬅️ **Negative (–)** → **Feature B → Feature A**
        """
    )

# Check if metadata is loaded
if 'md_for_analysis' in st.session_state and not st.session_state['md_for_analysis'].empty:
    
    chemprop_md = st.session_state['md_for_analysis']

    df = inside_levels(chemprop_md)
    mask = df.apply(lambda row: len(row['LEVELS']) == 0, axis=1)
    df = df[~mask]
    
    # Display the overview dataframe and checkboxes
    st.markdown("### Select Metadata Columns for ChemProp2 Calculation")

    st.markdown("""
                - **ChemProp2**: Select the metadata column that contains the timepoint information. This column must be present in the metadata for ChemProp2 calculation.
                - **Treatments**: Select the metadata column that contains the treatment conditions (e.g., control and treatment groups). This category is optional.
                - **Replicates**: Select the metadata column that contains the replicate information. This category is optional and will be used for plotting intensity trends.
                """)
    
    st.markdown("""
                <div style="color: #05577C; font-size: 20px;">
                Select <b>only one metadata attribute</b> per checkbox column.
                </div>""", 
                unsafe_allow_html=True)

    # Add a new column for checkboxes (initialized to False)
    df['ChemProp2 (Timepoint)'] = False
    df['Treatments'] = False
    df['Replicates'] = False

    # Display an editable DataFrame using st.data_editor
    edited_df = st.data_editor(df, num_rows="fixed", 
                               use_container_width=True, 
                               hide_index = True,
                               disabled=("ATTRIBUTES", "LEVELS", "COUNTS"),
                               column_order=("ChemProp2 (Timepoint)", "Treatments" , "Replicates", "ATTRIBUTES", "LEVELS", "COUNTS")
                               )

    # Collect the rows where checkboxes are selected (True)
    chemprop2_row = edited_df[edited_df['ChemProp2 (Timepoint)'] == True]['ATTRIBUTES'].tolist()
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
    st.session_state['chemprop_ft'] = final_ft
    st.session_state['chemprop_subset_md'] = subset_chemprop_md
    st.session_state['chemprop_md'] = time_md

    st.markdown("### ChemProp2 Metadata with Selected Columns")
    tab_info= [
            ("chemprop_subset_md", "Filtered Metadata"),
            ("chemprop_ft", "Filtered Feature Table"),
        ]
    
    display_tables_in_tabs(tab_info)

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

        st.session_state['chemprop_ft'] = subgroup
        st.session_state['chemprop_subset_md'] = subgroup_md
        st.session_state['chemprop_md'] = subgroup_md[chemprop2_row]

        filtered_tab_info= [
            ("chemprop_subset_md", "Filtered Metadata"),
            ("chemprop_ft", "Filtered Feature Table"),
        ]

        # Display the filtered data immediately
        with st.expander(f"Filtered ChemProp2 Tables"):
            display_tables_in_tabs(filtered_tab_info)

    else:
        st.warning("No group selected.")       


# Streamlit Interface
st.markdown('## ChemProp2 Scores')

dfs = get_chemprop2_inputs_from_state()
gnps_task_id = (st.session_state.get("gnps_task_id") or "").strip()

if dfs is None:
    st.info("Upload Edge Table + Feature Table + Metadata to run ChemProp2.")
    st.stop()

network_df, features_df, metadata_df = dfs

run_flag, mode = chemprop2_controls()
if not run_flag:
    st.stop()

if "gnps_task_id" in st.session_state and st.session_state["gnps_task_id"]:
    st.session_state["gnps_task_id"] = str(st.session_state["gnps_task_id"]).strip()

if "last_chemprop2_mode" not in st.session_state:
    st.session_state["last_chemprop2_mode"] = mode

if st.session_state["last_chemprop2_mode"] != mode:
    st.session_state["show_filters_chemprop2"] = False  # reset checkbox
    st.session_state["last_chemprop2_mode"] = mode
    st.rerun()

if mode == "Provided Edge Table":
        chemprop2_result_df = run_chemprop2_pipeline(network_df, features_df, metadata_df)
       
        if chemprop2_result_df is None or chemprop2_result_df.empty:
            st.warning("ChemProp2 returned no results after filtering.")
            st.stop()

        render_chemprop2_results_summary(network_df, chemprop2_result_df)
        render_chemprop2_download_buttons(chemprop2_result_df)

        st.markdown('### False Discovery Rate')

        if st.button("Apply FDR"):
             if 'ChemProp2_scores' in st.session_state:
                st.write('Select the positive and negative cutoffs for the ChemProp2 scores based on your FDR-curve')
                overall_fdr_table, fig_histogram, fig_fdr = calculate_fdr(
                    network_df, 
                    features_df, 
                    metadata_df, 
                    score_range=(-1, 1), 
                    bin_size=0.001
                    )
                 # Sample and blank selection interface
                c1, c2 = st.columns(2)
                c1.plotly_chart(fig_histogram)
                c2.plotly_chart(fig_fdr)

        render_chemprop2_scores_plot(chemprop2_result_df)

        render_chemprop2_filters_and_plots(
            scores_df=chemprop2_result_df,
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

    chemprop2_result_df = run_chemprop2_pipeline(cascade_df, features_df, metadata_df)
    if chemprop2_result_df is None or chemprop2_result_df.empty:
            st.warning("ChemProp2 returned no results after filtering.")
            st.stop()

    render_chemprop2_results_summary(cascade_df, chemprop2_result_df)
    render_chemprop2_download_buttons(chemprop2_result_df)

    st.markdown('### False Discovery Rate')

    if st.button("Apply FDR"):
         if 'ChemProp2_scores' in st.session_state:
            st.write('Select the positive and negative cutoffs for the ChemProp2 scores based on your FDR-curve')
            overall_fdr_table, fig_histogram, fig_fdr = calculate_fdr(
                cascade_df, 
                features_df, 
                metadata_df, 
                score_range=(-1, 1), 
                bin_size=0.001
                )
             # Sample and blank selection interface
            c1, c2 = st.columns(2)
            c1.plotly_chart(fig_histogram)
            c2.plotly_chart(fig_fdr)

    render_chemprop2_scores_plot(chemprop2_result_df)
    render_chemprop2_filters_and_plots(
        scores_df=chemprop2_result_df,
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

    chemprop2_result_df = run_chemprop2_pipeline(user_edge_df, features_df, metadata_df)

    if chemprop2_result_df is None or chemprop2_result_df.empty:
        st.warning("ChemProp2 returned no results for the selected edge (feature missing after filtering?).")
        st.stop()

    render_chemprop2_results_summary(user_edge_df, chemprop2_result_df)
    render_chemprop2_edge_detail_plots(
        filtered_df=chemprop2_result_df,
        all_scores_df=chemprop2_result_df,
        features_df=features_df,
        metadata_df=metadata_df
        )
    if (st.session_state.get("gnps_task_id") or "").strip():
        render_spectra_modifinder(chemprop2_result_df, st.session_state.get("an_gnps"))




        