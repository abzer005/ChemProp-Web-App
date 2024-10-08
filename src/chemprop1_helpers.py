import streamlit as st
from .common import *  # Importing common functionalities from the 'common' module
from src.cleanup import *
import pandas as pd
import re
import numpy as np
from scipy.stats import spearmanr, pearsonr

def subset_levels(metatable, feature_table):
    if metatable.empty or feature_table.empty:
        st.warning("Metadata or Feature table is empty. Please check the data before subsetting.")
        return metatable, feature_table

    # User selects the attributes to subset
    all_columns = metatable.columns.tolist()
    selected_attributes = st.multiselect("Select multiple attributes to subset:", all_columns, default=all_columns[0])

    for attribute in selected_attributes:
        levels_cdtn = metatable[attribute].unique()  # unique levels for the selected attribute

        # User decision to keep or exclude conditions
        c1,c2 = st.columns(2)
        read_cdtn = c1.radio(f"Do you want to keep or exclude conditions for '{attribute}'?", 
                             ('Keep', 'Exclude'), 
                             key = f'radio_{attribute}', )

        # User selects conditions to keep or exclude
        selected_conditions = c2.multiselect(f"Select conditions to {read_cdtn.lower()} for '{attribute}':", 
                                             levels_cdtn, 
                                             key = f'multiselect_{attribute}')

        if len(selected_conditions) > 0:
             # Filter the metatable based on the user's selection
            if read_cdtn == 'Keep':
                metatable = metatable[metatable[attribute].isin(selected_conditions)]
            elif read_cdtn == 'Exclude':
                metatable = metatable[~metatable[attribute].isin(selected_conditions)]

    # Check if subset_metadata index is in feature table columns
    if set(metatable.index).issubset(set(feature_table.columns)):
        # Subset the feature table based on the subset metadata index
        subset_feature_table = feature_table[metatable.index]

        with st.expander(f"Metadata after subsetting {metatable.shape}"):
            st.dataframe(metatable)
        
        with st.expander(f"Feature table after subsetting {subset_feature_table.shape}"):
            st.dataframe(subset_feature_table)
    
    else:
        # Handle the case where some metadata indices are not in feature table columns
        missing_columns = set(metatable.index) - set(feature_table.columns)
        st.warning(f"Some metadata indices are not present as columns in the feature table: {missing_columns}")
        subset_feature_table = feature_table.copy()

    # Return the subsetted metadata and feature table
    return metatable, subset_feature_table


# Function to check if an attribute contains numerical values
def contains_convertible_values(series):
    for value in series.unique():
        # Try to convert each unique value to a number
        numeric_part = re.sub(r"[^\d.]", "", str(value))
        numeric_value = float(numeric_part) if numeric_part else 0.1

         # If the numeric value is zero, change it to 0.1
        if numeric_value == 0:
            numeric_value = 0.1
        try:
            float(numeric_part)
            return True  # Return True if conversion is successful
        except ValueError:
            continue  # Skip if conversion fails
    return False  # Return False if no values could be converted


def convert_series_values(series):
    converted_series = series.copy()  # Copy the series to avoid modifying the original one
    for index, value in series.items():
        # Strip away non-numeric characters
        numeric_part = re.sub(r"[^\d.]", "", str(value))

        # Convert to float; if conversion fails or results in zero, set to 0.0
        try:
            numeric_value = float(numeric_part)
            if numeric_value == 0:
                numeric_value = 0.0
        except ValueError:
            numeric_value = 0.0

        # Update the value in the copied series
        converted_series.at[index] = numeric_value

    return converted_series

#Chemprop analysis:

def run_chemprop1_analysis(network_table, metadata, feature_table, chemprop_attribute):
    
    network_table['Chemprop_scores'] = np.nan  # Initialize with NaN

    # Iterate over each row in the network table
    for i in range(len(network_table)):

        cluster_id1 = network_table.iloc[i, 0]
        cluster_id2 = network_table.iloc[i, 1]

        if cluster_id1 in feature_table.index and cluster_id2 in feature_table.index:
            FeatureA = feature_table.loc[cluster_id1]
            FeatureB = feature_table.loc[cluster_id2]

            # Get timepoints
            timepoints = metadata[chemprop_attribute].unique()
            timepoint1, timepoint2 = timepoints[0], timepoints[1]

            # Subset the features based on timepoints
            A1 = FeatureA[metadata.index[metadata[chemprop_attribute] == timepoint1]]
            A2 = FeatureA[metadata.index[metadata[chemprop_attribute] == timepoint2]]
            B1 = FeatureB[metadata.index[metadata[chemprop_attribute] == timepoint1]]
            B2 = FeatureB[metadata.index[metadata[chemprop_attribute] == timepoint2]]

            k = 1.0e-10
            score = np.log(((A1 + k) / (B1 + k)) / ((A2 + k) / (B2 + k)))
            st.write(score)

            # Handling possible NaN values and ensuring a single score per pair
            score = np.nanmean(score) if not np.isnan(score_values).all() else np.nan

            network_table.at[i, 'Chemprop_scores'] = score

            network_table.at[i, 'Chemprop_scores'] = score
        else:
            network_table.at[i, 'Chemprop_scores'] = 'NA'

    return network_table


# Function to strip non-numeric characters
def strip_non_numeric(val):
    if isinstance(val, str):  # Only apply the regex to string values
        return re.sub(r'[^\d.]+', '', val)  # Remove everything except digits and dots
    return val  # If it's not a string (e.g., already numeric), return it as-is


def ChemProp2_scoring(desired_network_table, desired_feature_table, names_nodes, chemprop_metadata):
    # Initialize the score columns with NA values (using np.nan for missing values)
    ChemProp2 = np.full(len(desired_network_table), np.nan)
    ChemProp_spearman = np.full(len(desired_network_table), np.nan)
    ChemProp_log = np.full(len(desired_network_table), np.nan)
    ChemProp_sqrt = np.full(len(desired_network_table), np.nan)
    
    for i in range(len(desired_network_table)):
        try:
            # Attempt to find the indices of the feature IDs
            feature_ids = [names_nodes.index(desired_network_table.iloc[i, 0]), 
                           names_nodes.index(desired_network_table.iloc[i, 1])]
            
            # Subset feature data
            x = desired_feature_table.iloc[feature_ids, :]
            
            # Get metadata timepoints
            chemprop_sample_names = chemprop_metadata.index
            chemprop_timepoint = chemprop_metadata.iloc[:, 0]
            
            # Reorder x based on sample names in metadata
            reorder_id = [list(desired_feature_table.columns).index(name) for name in chemprop_sample_names]
            reordered_x = x.iloc[:, reorder_id].T
            reordered_x = pd.concat([chemprop_timepoint.reset_index(drop=True), reordered_x.reset_index(drop=True)], axis=1)
            
            # Pearson correlation
            corr_result = reordered_x.corr(method="pearson")
            ChemProp2[i] = (corr_result.iloc[0, 2] - corr_result.iloc[0, 1]) / 2
            
            # Spearman correlation
            corr_spearman = reordered_x.corr(method="spearman")
            ChemProp_spearman[i] = (corr_spearman.iloc[0, 2] - corr_spearman.iloc[0, 1]) / 2
            
            # Log transformation
            log_reorderedX = pd.concat([reordered_x.iloc[:, 0], np.log1p(reordered_x.iloc[:, 1:3])], axis=1)
            corr_log = log_reorderedX.corr(method="pearson")
            ChemProp_log[i] = (corr_log.iloc[0, 2] - corr_log.iloc[0, 1]) / 2
            
            # Square root transformation
            sqrt_reorderedX = pd.concat([reordered_x.iloc[:, 0], np.sqrt(reordered_x.iloc[:, 1:3])], axis=1)
            corr_sqrt = sqrt_reorderedX.corr(method="pearson")
            ChemProp_sqrt[i] = (corr_sqrt.iloc[0, 2] - corr_sqrt.iloc[0, 1]) / 2
        
        except ValueError:
            # If the feature ID is not found, leave the ChemProp scores as np.nan (NA in Pandas)
            ChemProp2[i] = np.nan
            ChemProp_spearman[i] = np.nan
            ChemProp_log[i] = np.nan
            ChemProp_sqrt[i] = np.nan

    # Combine results into a new DataFrame
    new_network_table = desired_network_table.copy()
    new_network_table['ChemProp2'] = ChemProp2
    new_network_table['ChemProp_spearman'] = ChemProp_spearman
    new_network_table['ChemProp_log'] = ChemProp_log
    new_network_table['ChemProp_sqrt'] = ChemProp_sqrt

    # Add absolute value columns
    abs_values = new_network_table[['ChemProp2', 'ChemProp_spearman', 'ChemProp_log', 'ChemProp_sqrt']].abs()
    abs_values.columns = [f"abs_{col}" for col in abs_values.columns]
    
    # Add sign of ChemProp2
    new_network_table['Sign_ChemProp2'] = np.sign(new_network_table['ChemProp2'])
    
    # Concatenate absolute values
    ChemProp2_file = pd.concat([new_network_table, abs_values], axis=1)

    return ChemProp2_file

