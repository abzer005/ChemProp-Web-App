# Import necessary libraries
import streamlit as st
from src.chemprop1_helpers import *
from src.chemprop2_helpers import *
import pandas as pd
import re
import plotly.express as px

page_setup()

st.markdown("## ChemProp2 Analysis for Time Series Data")

st.info("""
**To run ChemProp2:**
        
Ensure that your metadata table is **properly formatted**, including: 
  - A **filename** column matching the feature table
  - A valid time column with numeric values for the timepoints (>2 timepoints).
  - Please have the timepoint categories arranged in the desired analysis order (e.g., increasing time), as this order is used during processing and is required for the workflow to run correctly.

👉 Please refer to the **Data Preparation Essentials** section on the **Home** page for
details on required metadata structure and formatting.
""")

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

    st.info("""
    - **ChemProp2 Timepoint (required):**  
    In the following checkbox table, select the metadata column that contains the **time information**.  
    This column contains the time information (e.g., 0, 12, 24 h) and is required to run ChemProp2.

    - **Treatments (optional):**  
    Select the metadata column that defines treatment groups (e.g., control vs treatment).
    This can be used for filtering and performing group-wise analyses.

    - **Replicates (optional):**  
    Select the metadata column that defines biological or technical replicates.  
    This information is used for **visualizing intensity trends**.
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

    # Enforce exactly one timepoint column for ChemProp2
    timepoint_column = chemprop2_row[0] if chemprop2_row else None

    if not timepoint_column:
        st.warning("Please select exactly one metadata column for ChemProp2 (Timepoint).")
        st.stop()

    # Keep ALL metadata columns for later filtering (do this once)
    st.session_state["chemprop2_md_full"] = chemprop_md.copy()

    # Store the user's "role assignments" (for later plotting/labels; optional)
    st.session_state["chemprop2_time_col"] = timepoint_column
    st.session_state["chemprop2_treatment_cols"] = treatment_row
    st.session_state["chemprop2_replicate_cols"] = rep_row

    # Reset ChemProp2 filters if timepoint column changed
    prev_tp = st.session_state.get("_chemprop2_prev_time_col")
    if prev_tp != timepoint_column:
        st.session_state.pop("chemprop2_extra_filters", None)
        st.session_state["_chemprop2_prev_time_col"] = timepoint_column

    # Show the selected attributes for Treatments
    if chemprop2_row:
        st.success(f"The column selected for calculating ChemProp2: {', '.join(chemprop2_row)}")
    else:
        st.warning("No column selected for calculating ChemProp2 scores.")

    # Show the selected attributes for Treatments
    if treatment_row:
        st.success(f"The column selected for Treatment conditions: {', '.join(treatment_row)}")
    else:
        st.warning("No column selected for Treatments. This is optional and you can proceed without it.")

    # Show the selected attributes for Replicates
    if rep_row:
        st.success(f"The column selected for Replicate conditions: {', '.join(rep_row)}")
    else:
        st.warning("No column selected for Replicates. This is optional and you can proceed without it.")
    
    # Apply the function to the subset DataFrame
    time_md = (
        chemprop_md[[timepoint_column]]
        .applymap(strip_non_numeric)
        .apply(pd.to_numeric, errors="coerce")
    ) if timepoint_column else pd.DataFrame(index=chemprop_md.index)

    subset_chemprop_md = chemprop_md.copy()
    subset_chemprop_md[timepoint_column] = time_md[timepoint_column]

    st.markdown("### ChemProp2 Metadata with Selected Columns")
        
    final_ft = st.session_state['ft_for_analysis'].copy()
    st.session_state['chemprop2_ft'] = final_ft
    st.session_state['chemprop2_subset_md'] = subset_chemprop_md
    st.session_state["chemprop2_md"] = subset_chemprop_md[[timepoint_column]] if timepoint_column else pd.DataFrame()

    tab_info= [
            ("chemprop2_subset_md", "Filtered Metadata"),
            ("chemprop2_ft", "Filtered Feature Table"),
        ]
    
    display_tables_in_tabs(tab_info)

else:
    st.warning("Metadata not loaded in the session state.")


st.markdown("### Filter the data")

with st.expander("ℹ️ Why would I want to filter the data? (optional)"):
    st.markdown("""
Filtering is **optional**, but can be helpful depending on your experimental design and analysis goals.
Using this section, you can optionally subset by **any metadata column** (up to 10 filters).  

You may want to filter your data to:
- Restrict the analysis to **specific timepoints** of interest
- Run ChemProp **separately for different subgroups**
- Focus on a subset of samples or replicates
                
For example, when both **treatment** and **control** groups are present, it is recommended to run
ChemProp2 **independently for each group** rather than analyzing them together.
""")
    
do_subset = st.checkbox("Do you want to subset the data further?", False)
if do_subset:
    st.markdown("Select samples based on metadata filters.")

    MAX_EXTRA_FILTERS = 10
    ss_filters_key = "chemprop2_extra_filters"

    if ss_filters_key not in st.session_state:
        st.session_state[ss_filters_key] = []  # [{"column": None, "values": []}, ...]

    # Start from metadata that contains ALL columns
    md_base = st.session_state.get("chemprop1_md_full", subset_chemprop_md).copy()

    if "chemprop2_md" in st.session_state and isinstance(st.session_state["chemprop2_md"], pd.DataFrame):
        if timepoint_column in st.session_state["chemprop2_md"].columns:
            md_base[timepoint_column] = st.session_state["chemprop2_md"][timepoint_column]

    # Feature table aligned to metadata samples
    ft_base = st.session_state["ft_for_analysis"].copy()
    ft_base = ft_base.loc[:, md_base.index]  # ensure alignment

    st.write("#### Additional Filters (optional)")

    col_add, col_clear = st.columns([1, 1])
    with col_add:
        if st.button("➕ Add filter", disabled=(len(st.session_state[ss_filters_key]) >= MAX_EXTRA_FILTERS)):
            st.session_state[ss_filters_key].append({"column": None, "values": []})

    with col_clear:
        if st.button("🧹 Clear all filters", disabled=(len(st.session_state[ss_filters_key]) == 0)):
            st.session_state[ss_filters_key] = []

    available_columns = list(md_base.columns)

    for i, f in enumerate(st.session_state[ss_filters_key]):
        with st.container():
            left, mid, right = st.columns([3, 5, 1])

            selected_col = left.selectbox(
                f"Filter {i+1} – Column",
                options=["(choose)"] + available_columns,
                index=0 if not f.get("column") else (available_columns.index(f["column"]) + 1),
                key=f"chemprop2_filter_col_{i}",
            )
            selected_col = None if selected_col == "(choose)" else selected_col
            st.session_state[ss_filters_key][i]["column"] = selected_col

            if selected_col:
                val_options = sorted(md_base[selected_col].dropna().astype(str).unique())
            else:
                val_options = []

            selected_vals = mid.multiselect(
                f"Filter {i+1} – Values",
                options=val_options,
                default=[v for v in (f.get("values") or []) if str(v) in val_options],
                key=f"chemprop2_filter_vals_{i}",
            )
            st.session_state[ss_filters_key][i]["values"] = selected_vals

            if right.button("✖", key=f"chemprop2_filter_remove_{i}"):
                st.session_state[ss_filters_key].pop(i)
                st.rerun()

    # Apply filters sequentially
    final_md, final_ft = apply_md_filters(md_base, ft_base, st.session_state[ss_filters_key])

    # Save to session state for later use
    st.session_state["chemprop2_ft"] = final_ft
    st.session_state["chemprop2_subset_md"] = final_md

    timepoint_column = st.session_state.get("chemprop2_time_col")  
    st.session_state["chemprop2_md"] = final_md[[timepoint_column]] if timepoint_column else pd.DataFrame()
    
    filtered_tab_info = [
        ("chemprop2_subset_md", "Filtered Metadata"),
        ("chemprop2_ft", "Filtered Feature Table"),
    ]

    with st.expander("Filtered ChemProp2 Tables", expanded=False):
        display_tables_in_tabs(filtered_tab_info)

    st.caption(f"Samples remaining: {final_md.shape[0]}")
    
    if final_md.shape[0] == 0:
        st.warning("No samples left after filtering — loosen one of the filters.")   


# Streamlit Interface
st.markdown('## ChemProp2 Scores')
#st.dataframe(st.session_state["chemprop2_md"].head(), use_container_width=True)

dfs = get_chemprop2_inputs_from_state()
gnps_task_id = (st.session_state.get("gnps_task_id") or "").strip()

if dfs is None:
    st.info("Upload Edge Table + Feature Table + Metadata to run ChemProp2.")
    st.stop()

network_df, features_df, metadata_df = dfs

time_series = st.session_state.get("chemprop2_md")

if isinstance(time_series, pd.DataFrame) and not time_series.empty:
    time_col = time_series.columns[0]
    if time_col in metadata_df.columns:
        metadata_df = sort_metadata_by_selected_time(metadata_df, time_col)

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

                with st.expander("ℹ️ How to use the FDR plots to choose ChemProp2 cutoffs"):
                    st.info("""
                ChemProp2 provides an **FDR-based guide** to help you choose reliable **positive and negative score cutoffs**.

                - **The histogram on the left** compares ChemProp2 scores using the **original feature table (target)** with scores from a
                **randomized decoy feature table**. This shows which score ranges are likely to occur by chance.  
                Some overlap between the two is expected and **not a problem**. What matters most is the
                **tail of the distribution**, where target scores clearly exceed decoy scores.

                - **The FDR curve on the right** is calculated from cumulative FDR values across score bins. 
                Use the FDR curve to select score thresholds that match your confidence level.
                For example, a **1% FDR** may suggest keeping only scores above **±0.6**, while a **5% FDR**
                may allow more relaxed cutoffs (e.g., **±0.4**).

                Once selected, these cutoffs can be applied to filter results from your downloaded ChemProp2 score table.
                """)

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
            
            with st.expander("ℹ️ How to use the FDR plots to choose ChemProp2 cutoffs"):
                    st.info("""
                ChemProp2 provides an **FDR-based guide** to help you choose reliable **positive and negative score cutoffs**.

                - **The histogram on the left** compares ChemProp2 scores using the **original feature table (target)** with scores from a
                **randomized decoy feature table**. This shows which score ranges are likely to occur by chance.  
                Some overlap between the two is expected and **not a problem**. What matters most is the
                **tail of the distribution**, where target scores clearly exceed decoy scores.

                - **The FDR curve on the right** is calculated from cumulative FDR values across score bins. 
                Use the FDR curve to select score thresholds that match your confidence level.
                For example, a **1% FDR** may suggest keeping only scores above **±0.6**, while a **5% FDR**
                may allow more relaxed cutoffs (e.g., **±0.4**).

                Once selected, these cutoffs can be applied to filter results from your downloaded ChemProp2 score table.
                """)
                    
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




        