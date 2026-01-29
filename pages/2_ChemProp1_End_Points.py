# Import necessary libraries
import streamlit as st
from src.chemprop1_helpers import *
import pandas as pd
import streamlit_agraph

page_setup()

st.markdown("### ChemProp1 Analysis")

st.warning(
    """
- ChemProp1 is designed for experiments with **exactly 2 sequential timepoints**. 
- If your dataset contains more than 2 timepoints, please use the **Filter the Data** 
section below to subset your data to **2 timepoints only**, or proceed to 
**ChemProp2** for multi-timepoint analysis.

""")

st.info("""
**To run ChemProp1:**
        
Ensure that your metadata table is **properly formatted**, including: 
  - A **filename** column matching the feature table
  - A valid time column with numeric values for the two time points.
  - Please have the timepoint categories arranged in the desired analysis order (e.g., increasing time), as this order is used during processing and is required for the workflow to run correctly.

👉 Please refer to the **Data Preparation Essentials** section on the **Home** page for
details on required metadata structure and formatting.
""")

# Check if metadata is loaded
if 'md_for_analysis' in st.session_state and not st.session_state['md_for_analysis'].empty:
    
    chemprop_md = st.session_state['md_for_analysis']

    df = inside_levels(chemprop_md)
    mask = df.apply(lambda row: len(row['LEVELS']) == 0, axis=1)
    df = df[~mask]
    
    # Display the overview dataframe and checkboxes
    st.markdown("### Select Metadata Columns for ChemProp1 Calculation")

    st.info("""
- **ChemProp1 Timepoint (required):**  
  In the following checkbox table, select the metadata column that contains the **time information**.  
  This column **must contain exactly two distinct values** and is required to run ChemProp1.

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

    # Enforce exactly one timepoint column for ChemProp1
    timepoint_column = chemprop_row[0] if chemprop_row else None

    # Keep ALL metadata columns for later filtering (do this once)
    st.session_state["chemprop1_md_full"] = chemprop_md.copy()

    # Store the user's "role assignments" (for later plotting/labels; optional)
    st.session_state["chemprop1_time_col"] = timepoint_column
    st.session_state["chemprop1_treatment_cols"] = treatment_row
    st.session_state["chemprop1_replicate_cols"] = rep_row


    # Show the selected attributes for Treatments
    if chemprop_row:
         # Extract unique values from the selected attributes in `chemprop_md` dataframe
        if timepoint_column:
            unique_count = chemprop_md[timepoint_column].nunique(dropna=True)

            if unique_count < 3:
                st.success(
                    f"You have {unique_count} unique timepoints in '{timepoint_column}'. "
                    "ChemProp1 is appropriate for this data.",
                    icon="✅"
                )
            else:
                st.error(
                    "You cannot run ChemProp1 for attributes with more than 2 unique timepoints. "
                    "Consider running ChemProp2, or filter the data further in the 'Filter the Data' section below to have only 2 timepoints.",
                    icon="🚨"
                )

            st.success(f"The column selected for calculating ChemProp1: {timepoint_column}")
        else:
            st.warning("No column selected for calculating ChemProp1 scores.")

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
    
    # Subset the DataFrames based on the user's selection
    treatment_md = chemprop_md[treatment_row] if treatment_row else pd.DataFrame()
    rep_md = chemprop_md[rep_row] if rep_row else pd.DataFrame()
    
    time_md = (
        chemprop_md[[timepoint_column]]
        .applymap(strip_non_numeric)
        .apply(pd.to_numeric, errors="coerce")
    ) if timepoint_column else pd.DataFrame(index=chemprop_md.index)

    # Keep ALL metadata columns for later filtering
    st.session_state["chemprop1_md_full"] = chemprop_md.copy()

    subset_chemprop_md = chemprop_md.copy()
    if timepoint_column:
        subset_chemprop_md[timepoint_column] = time_md[timepoint_column]
        
    st.markdown("### ChemProp1 Metadata with Selected Columns")

    final_ft = st.session_state['ft_for_analysis'].copy()
    st.session_state['chemprop1_ft'] = final_ft
    st.session_state['chemprop1_subset_md'] = subset_chemprop_md
    st.session_state["chemprop1_md"] = subset_chemprop_md[[timepoint_column]] if timepoint_column else pd.DataFrame()

    tab_info= [
        ("chemprop1_subset_md", "Filtered Metadata"),
        ("chemprop1_ft", "Filtered Feature Table"),
    ]

    display_tables_in_tabs(tab_info)

else:
    st.warning("Metadata not loaded in the session state.")

################################## FILTER THE DATA ##########################################################
st.markdown("### Filter the Data")

with st.expander("ℹ️ Why would I want to filter the data? (optional)"):
    st.markdown("""
Filtering is **optional**, but can be helpful depending on your experimental design and analysis goals.

You may want to filter your data to:
- Restrict the analysis to **specific timepoints** (e.g., required for ChemProp1, only 2 timepoints)
- Run ChemProp **separately for different subgroups**
- Focus on a subset of samples or replicates

This section provides:  
1️⃣ **Initial Timepoint Filter** – select **exactly two timepoints** (required for ChemProp1).  
2️⃣ **Additional Filters** – optionally subset by **any metadata column** (up to 10 filters).  
   For example, if your dataset includes **treatment** and **control** groups, it is recommended to run 
                ChemProp **independently for each group** rather than combining them.
""")

do_subset = st.checkbox("Do you want to subset the data further?", False)

if do_subset:
    st.markdown("Select samples based on their conditions.")

    # --- Determine the (single) timepoint column ---
    timepoint_column = st.session_state.get("chemprop1_time_col")
    if not timepoint_column:
        # fallback if session_state not set yet
        timepoint_column = chemprop_row[0] if isinstance(chemprop_row, list) and chemprop_row else chemprop_row

    if not timepoint_column:
        st.warning("Please select a timepoint column above before filtering.")
        st.stop()

    # Base metadata for filtering
    md_base = st.session_state.get("chemprop1_md_full", subset_chemprop_md).copy()

    # Make sure the cleaned numeric timepoint column is used for ChemProp1 logic
    if "chemprop1_md" in st.session_state and isinstance(st.session_state["chemprop1_md"], pd.DataFrame):
        if timepoint_column in st.session_state["chemprop1_md"].columns:
            md_base[timepoint_column] = st.session_state["chemprop1_md"][timepoint_column]

    # Feature table aligned to metadata samples
    ft_base = st.session_state["ft_for_analysis"].copy()
    ft_base = ft_base.loc[:, md_base.index]  # ensure alignment

    # -------------------------------
    # 1) Initial Timepoint Filter (exactly 2)
    # -------------------------------
    with st.container():
        st.write("#### Initial Timepoint Filter")
        c1, c2 = st.columns(2)

        c1.selectbox(
            "Select the column for categorization (Timepoints)",
            options=[timepoint_column],
            index=0,
            disabled=True,
        )

        subset_group = c2.multiselect(
            "Select exactly TWO timepoints for filtering",
            options=sorted(md_base[timepoint_column].dropna().astype(str).unique()),
            max_selections=2,
        )

    if len(subset_group) != 2:
        st.info("Choose exactly two timepoints to continue.")
        st.stop()

    subset_group = list(map(str, subset_group))
    tp_mask = md_base[timepoint_column].astype(str).isin(subset_group)

    md_after_tp = md_base.loc[tp_mask].copy()
    ft_after_tp = ft_base.loc[:, md_after_tp.index].copy()

    # -------------------------------
    # 2) Additional Filters (any column), up to 10
    # -------------------------------
    st.write("#### Additional Filters (optional)")

    MAX_EXTRA_FILTERS = 10
    if "chemprop1_extra_filters" not in st.session_state:
        st.session_state["chemprop1_extra_filters"] = []  # [{"column": None, "values": []}, ...]

    col_add, col_clear = st.columns([1, 1])

    with col_add:
        if st.button(
            "➕ Add filter",
            disabled=(len(st.session_state["chemprop1_extra_filters"]) >= MAX_EXTRA_FILTERS),
        ):
            st.session_state["chemprop1_extra_filters"].append({"column": None, "values": []})

    with col_clear:
        if st.button("🧹 Clear all filters", disabled=(len(st.session_state["chemprop1_extra_filters"]) == 0)):
            st.session_state["chemprop1_extra_filters"] = []

    available_columns = list(md_after_tp.columns)
    # Optional: don't offer the timepoint column again since it's already handled
    # available_columns = [c for c in md_after_tp.columns if c != timepoint_column]

    for i, f in enumerate(st.session_state["chemprop1_extra_filters"]):
        with st.container():
            left, mid, right = st.columns([3, 5, 1])

            selected_col = left.selectbox(
                f"Filter {i+1} – Column",
                options=["(choose)"] + available_columns,
                index=0 if not f.get("column") else (available_columns.index(f["column"]) + 1),
                key=f"chemprop1_filter_col_{i}",
            )
            selected_col = None if selected_col == "(choose)" else selected_col
            st.session_state["chemprop1_extra_filters"][i]["column"] = selected_col

            if selected_col:
                val_options = sorted(md_after_tp[selected_col].dropna().astype(str).unique())
            else:
                val_options = []

            selected_vals = mid.multiselect(
                f"Filter {i+1} – Values",
                options=val_options,
                default=[v for v in (f.get("values") or []) if str(v) in val_options],
                key=f"chemprop1_filter_vals_{i}",
            )
            st.session_state["chemprop1_extra_filters"][i]["values"] = selected_vals

            if right.button("✖", key=f"chemprop1_filter_remove_{i}"):
                st.session_state["chemprop1_extra_filters"].pop(i)
                st.rerun()

    # Apply extra filters sequentially
    final_md, final_ft = apply_md_filters(md_after_tp, ft_after_tp, st.session_state["chemprop1_extra_filters"])

    # -------------------------------
    # Save results + show tables
    # -------------------------------
    st.session_state["chemprop1_ft"] = final_ft
    st.session_state["chemprop1_subset_md"] = final_md
    st.session_state["chemprop1_md"] = final_md[[timepoint_column]]

    with st.expander("Filtered Tables", expanded=False):
        display_tables_in_tabs([
            ("chemprop1_subset_md", "Filtered Metadata"),
            ("chemprop1_ft", "Filtered Feature Table"),
        ])

    st.caption(f"Samples remaining: {final_md.shape[0]}")
    if final_md.shape[0] == 0:
        st.warning("No samples left after filtering — loosen one of the filters.")


# Streamlit Interface
st.markdown('## Calculate ChemProp1 Scores')

dfs = get_inputs_from_state()

if dfs is None:
    st.info("Upload Edge Table + Feature Table + Metadata to run ChemProp1.")
    st.stop()

network_df, features_df, metadata_df = dfs

time_series = st.session_state.get("chemprop1_md")

if time_series is not None:
    time_col = time_series.columns[0]

    metadata_df = sort_metadata_by_selected_time(metadata_df, time_col)

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

