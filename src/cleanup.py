import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from scipy.spatial import distance
from sklearn.preprocessing import StandardScaler

# currently conflicting dependencies (requires old pandas 1.2.4)

@st.cache_data
def clean_up_md(md):
    md = (
        md.copy()
    )  # storing the files under different names to preserve the original files
    # remove the (front & tail) spaces, if any present, from the rownames of md
    md = md.dropna(how="all")
    md.index = [name.strip() for name in md.index]
    # for each col in md
    # 1) removing the spaces (if any)
    # 2) replace the spaces (in the middle) to underscore
    # 3) converting them all to UPPERCASE
    for col in md.columns:
        if md[col].dtype == str:
            md[col] = [item.strip().replace(" ", "_").upper() for item in md[col]]
    md.index = [i.replace(".mzXML", "").replace(".mzML", "").replace(" Peak area", "") for i in md.index]
    return md


@st.cache_data
def clean_up_ft(ft):
    ft = (
        ft.copy()
    )  # storing the files under different names to preserve the original files
    
    ft = ft.dropna(how="all")

    # drop all columns that are not mzML or mzXML file names
    ft.drop(columns=[col for col in ft.columns if not (".mzML" in col or ".mzXML" in col)], inplace=True)
    
    # remove " Peak area" from column names, contained after mzmine pre-processing
    ft.rename(
        columns={col: col.replace(" Peak area", "").replace(".mzXML", "").replace(".mzML", "").strip() for col in ft.columns},
        inplace=True,
    )
    return ft


@st.cache_data
def check_columns(md, ft):
    if sorted(ft.columns) != sorted(md.index):
        st.warning("Not all files are present in both meta data & feature table.")
        
        # Find and remove columns in 'ft' that are not in the index of 'md'
        ft_cols_not_in_md = [col for col in ft.columns if col not in md.index]
        if ft_cols_not_in_md:
            st.warning(
                f"These {len(ft_cols_not_in_md)} columns of feature table are not present in metadata table and will be removed:\n{', '.join(ft_cols_not_in_md)}"
            )
            ft = ft.drop(columns=ft_cols_not_in_md)
        
        # Find and remove rows in 'md' that are not in the columns of 'ft'
        md_rows_not_in_ft = [row for row in md.index if row not in ft.columns]
        if md_rows_not_in_ft:
            st.warning(
                f"These {len(md_rows_not_in_ft)} rows of metadata table are not present in feature table and will be removed:\n{', '.join(md_rows_not_in_ft)}"
            )
            md = md.drop(md_rows_not_in_ft)
    return md, ft


@st.cache_data
def inside_levels(df):
    df = pd.DataFrame(
        {
            "ATTRIBUTES": df.columns,
            "LEVELS": [set(df[col].dropna().astype(str).to_list()) for col in df],
            "COUNTS": [df[col].value_counts().to_list() for col in df],
        }
    )
    return df


@st.cache_data
def get_cutoff_LOD(df):
    # get the minimal value that is not zero (lowest measured intensity)
    return round(df.replace(0, np.nan).min(numeric_only=True).min())


@st.cache_data
def remove_blank_features(blanks, samples, cutoff):
    # Getting mean for every feature in blank and Samples
    avg_blank = blanks.mean(
        axis=1, skipna=False
    )  # set skipna = False do not exclude NA/null values when computing the result.
    avg_samples = samples.mean(axis=1, skipna=False)

    # Getting the ratio of blank vs samples
    ratio_blank_samples = (avg_blank + 1) / (avg_samples + 1)

    # Create an array with boolean values: True (is a real feature, ratio<cutoff) / False (is a blank, background, noise feature, ratio>cutoff)
    is_real_feature = ratio_blank_samples < cutoff

    # Calculating the number of background features and features present (sum(bg_bin) equals number of features to be removed)
    n_background = len(samples) - sum(is_real_feature)
    n_real_features = sum(is_real_feature)

    blank_removal = samples[is_real_feature.values]

    return blank_removal, n_background, n_real_features


@st.cache_data
def impute_missing_values(df, cutoff_LOD):
    # impute missing values (0) with a random value between one and lowest intensity (cutoff_LOD)
    if cutoff_LOD > 1:
        return df.apply(
            lambda x: [np.random.randint(1, cutoff_LOD) if v == 0 else v for v in x]
        )

@st.cache_data
def normalization(feature_df, meta_data_df):
    # Transpose feature table for processing
    feature_df = feature_df.T

    # Remove metadata rows that are not in the feature table
    md_rows_not_in_samples = [row for row in meta_data_df.index if row not in feature_df.index]
    md_samples = meta_data_df.drop(md_rows_not_in_samples)

    # Put feature table and metadata in the same order
    feature_df.sort_index(inplace=True)
    md_samples.sort_index(inplace=True)

    # Ensure sample names are the same between feature and metadata tables
    try:
        if not md_samples.index.equals(feature_df.index):
            st.warning("Sample names in feature and metadata table are NOT the same!")
    except ValueError as e:
        st.warning(f"Sample names cannot be compared. Error: {e}")
        return None, None

    # Normalize feature data (handling division by zero)
    normalized = feature_df.apply(lambda x: x / np.sum(x) if np.sum(x) != 0 else 0, axis=1)
    
    # Return a copy to avoid caching issues with mutable objects
    return md_samples.copy(), normalized.copy()

# can not hash pcoa
def permanova_pcoa(scaled, distance_matrix, attribute):
    # Scale the data
    #scaler = StandardScaler()
    #scaled = scaler.fit_transform(unscaled_data)

    # Create the distance matrix from the original data
    distance_matrix = skbio.stats.distance.DistanceMatrix(
        distance.squareform(distance.pdist(scaled.values, distance_matrix))
    )

    # After calculating the distance matrix
    #st.write(f"Distance matrix shape: {distance_matrix.shape}")

    # perform PERMANOVA test
    permanova = skbio.stats.distance.permanova(distance_matrix, attribute)
    permanova["R2"] = 1 - 1 / (
        1
        + permanova["test statistic"]
        * permanova["number of groups"]
        / (permanova["sample size"] - permanova["number of groups"] - 1)
    )
    # perfom PCoA
    pcoa = skbio.stats.ordination.pcoa(distance_matrix)

    return permanova, pcoa


# can not hash pcoa
def get_pcoa_scatter_plot(pcoa, md_samples, attribute):
    df = pcoa.samples[["PC1", "PC2"]]
    df = df.set_index(md_samples.index)
    df = pd.merge(
        df[["PC1", "PC2"]],
        md_samples[attribute].apply(str),
        left_index=True,
        right_index=True,
    )

    title = f"PRINCIPLE COORDINATE ANALYSIS"
    fig = px.scatter(
        df,
        x="PC1",
        y="PC2",
        template="plotly_white",
        width=600,
        height=400,
        color=attribute,
        hover_name=df.index,
    )

    fig.update_layout(
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={"text": title, "font_color": "#3E3D53"},
        xaxis_title=f"PC1 {round(pcoa.proportion_explained[0]*100, 1)}%",
        yaxis_title=f"PC2 {round(pcoa.proportion_explained[1]*100, 1)}%",
    )
    return fig

# can not hash pcoa
def get_pcoa_variance_plot(pcoa):
    # To get a scree plot showing the variance of each PC in percentage:
    percent_variance = np.round(pcoa.proportion_explained * 100, decimals=2)

    fig = px.bar(
        x=[f"PC{x}" for x in range(1, len(pcoa.proportion_explained) + 1)],
        y=percent_variance,
        template="plotly_white",
        width=500,
        height=400,
    )
    fig.update_traces(marker_color="#696880", width=0.5)
    fig.update_layout(
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={"text": "PCoA - VARIANCE", "x": 0.5, "font_color": "#3E3D53"},
        xaxis_title="principal component",
        yaxis_title="variance (%)",
    )
    return fig
