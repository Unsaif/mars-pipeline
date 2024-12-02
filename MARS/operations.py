import pandas as pd
import json
import os
import numpy as np
import re
import logging

logger = logging.getLogger('main.operations')

def remove_clades_from_taxaNames(merged_df, taxaSplit='; '):
    """ 
    Removes clade extensions from taxonomic names at any taxonomic level (if present)
    and sums counts of all clades of each taxa together.

    Args:
        merged_df (pd.DataFrame): The input DataFrame with taxonomic groups in the index.
        taxaSplit (string):       Seperator by which input taxonomic levels are seperated.
    
    Returns:
        grouped_df (pd.DataFrame): The input DataFrame with taxonomic groups in the index, without clade seperation anymore.
    """
    logger.info("Removing clade extensions from taxa names.")
    
    merged_df = merged_df.reset_index()
    taxa = merged_df['Taxon']
    
    # Abbreviations for all taxonomic levels
    taxUnitAbbreviations = ['p__', 'c__', 'o__', 'f__', 'g__', 's__']
    # Specifies the column in which taxonomic level is present in split name (increases by 1 with each for-loop iteration)
    columnTracker = 1

    for taxLevel in taxUnitAbbreviations:
        # Filter taxa that contain taxLevel & split the strings by ';'
        taxa_wTaxLevel = taxa[taxa.str.contains(taxLevel)]
        taxaSep = taxa_wTaxLevel.str.split(taxaSplit, expand=True)
    
        # Extract taxa names from the taxLevel column
        taxName = taxaSep[columnTracker].astype(str)
        
        # Filter taxName with occourence of clade extensions to output in logger
        clade_pattern = r'(_clade)?_[a-zA-Z]$|(_clade)?_[a-zA-Z]\s'        
        taxNames_wClades = taxName.str.match(clade_pattern)
        matched_taxNames = taxName.loc[taxNames_wClades].tolist()        
        if matched_taxNames:
            logger_output = ', '.join(matched_taxNames)
            logger.info(f"List of {taxLevel} taxa for which clade extension was found & removed: {logger_output}.")
        else:
            logger.info(f"For {taxLevel} taxa no clade extension was found & removed.")

        # Clean taxa names using regex
        taxName = taxName.str.replace(r'(_clade)?_[a-zA-Z]$', '', regex=True)
        taxName = taxName.str.replace(r'(_clade)?_[a-zA-Z]\s', '', regex=True)
                
        # Update full taxa name (inlcuding all taxLevels) with cleaned taxa
        if taxLevel == 's__':
            taxaUpdate = pd.concat([taxaSep.iloc[:, :columnTracker], taxName], axis=1)
        else:
            taxaUpdate = pd.concat([taxaSep.iloc[:, :columnTracker], taxName, taxaSep.iloc[:, columnTracker+1:]], axis=1)
        taxaUpdate = taxaUpdate.apply(lambda x: taxaSplit.join(x.astype(str)), axis=1)
    
        # Replace original taxa with updated taxa names & update the merged_df
        taxa[taxa.str.contains(taxLevel)] = taxaUpdate.values

        # Increase column to match corresponding taxonomic level for next for-loop iteration
        columnTracker += 1
    
    # Replace the clade-containing taxa names by the updated taxa names without clades
    merged_df['Taxon'] = taxa

    # Group by 'Taxon' and sum, remove 'GroupCount' column if it exists &
    # rename columns to remove 'sum_' prefix
    grouped_df = merged_df.groupby('Taxon').sum().reset_index()
    if 'GroupCount' in grouped_df.columns:
        grouped_df.drop(columns='GroupCount', inplace=True)
    grouped_df.columns = [col.replace('sum_', '') for col in grouped_df.columns]

    return grouped_df


def split_taxonomic_groups(merged_df, flagLoneSpecies=False, taxaSplit='; '):
    """
    Split the taxonomic groups in the index of the input DataFrame and create separate DataFrames for each taxonomic level.

    Args:
        merged_df (pd.DataFrame): The input DataFrame with taxonomic groups in the index, without clade seperation anymore.

    Returns:
        dict: A dictionary with keys as taxonomic levels and values as the corresponding DataFrames.
    """
    logger.info("Split taxonomic levels in seperate dataframes, each per level.")

    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

    # Replace all level indicators in the 'Taxon' column
    # Added drop=True argument as indexing of merged_df is reset within 
    # remove_clades_from_taxaNames function - JW
    merged_df = merged_df.reset_index(drop=True)
    
    merged_df['Taxon'] = merged_df['Taxon'].replace(".__", "", regex=True)

    # Reset the index and split the index column into separate columns for each taxonomic level
    taxonomic_levels_df = merged_df 
    taxonomic_split_df = taxonomic_levels_df['Taxon'].str.split(taxaSplit, expand=True)
    taxonomic_split_df = taxonomic_split_df.fillna('') # deals with cases of empty string instead of np.nan
    taxonomic_split_df.columns = levels

    # Concatenate genus and species names if both are present, otherwise leave species column unchanged
    if flagLoneSpecies:
        taxonomic_split_df['Species'] = taxonomic_split_df.apply(
            lambda row: row['Genus'] + '_' + row['Species'] if row['Species'] != '' else row['Species'],
            axis=1
    )

    # Concatenate the taxonomic_split_df and the abundance data from taxonomic_levels_df
    taxonomic_levels_df = pd.concat([taxonomic_split_df, taxonomic_levels_df.iloc[:, 1:]], axis=1)

    # Initialize a dictionary to store DataFrames for each taxonomic level
    taxonomic_dfs = {}

    # Iterate through the taxonomic levels and create a DataFrame for each level
    for level in levels:
        level_df = taxonomic_levels_df[[level] + list(taxonomic_levels_df.columns[len(levels):])]
        level_df = level_df.rename(columns={level: 'Taxon'})

        # Set the 'Taxon' column as the index and remove rows with an empty string in the 'Taxon' column
        level_df = level_df.set_index('Taxon')
        level_df = level_df.loc[level_df.index != '']

        # Add the DataFrame to the dictionary
        taxonomic_dfs[level] = level_df

    return taxonomic_dfs


def filter_samples_low_read_counts(taxonomic_dataframes, sample_read_counts_cutoff=1):
    """
    Filter all taxonomic abundance dataframes (containing absolute read counts) 
    for samples which contain less total read counts than a specified cutoff & exclude them
    from the output dataframes. This ensures there are no samples with read counts = 0, which
    could lead to downstream errors in the pipeline, as well as that the sequencing depth
    is high enough in each sample that saturation for OTU detection is reached.

    Args:
        taxonomic_dataframes (dict):        A dictionary with keys as taxonomic levels and values as the corresponding DataFrames.
        sample_read_counts_cutoff (int):    A cutoff for minimal read counts in a sample to be included in downstream analysis.
                                            Defaults to 1.

    Returns:
        filtered_taxonomic_dataframes (dict): A dictionary with keys as taxonomic levels and values as the corresponding DataFrames,
                                            only with samples containing more reads than the cutoff.
    """
    logger.info('Filtering taxonomic dataframes for samples with read counts below/equal to cutoff (currently 1 read).')

    filtered_taxonomic_dataframes = {}

    for level, df in taxonomic_dataframes.items():
        read_counts = df.sum()
        
        # Subset the level dataframe only including those samples with read count higher than cutoff
        samples_higher_than_cutoff = [i for i, v in enumerate(read_counts) if v > sample_read_counts_cutoff]
        level_filtered_taxonomic_dataframe = df.iloc[:, samples_higher_than_cutoff]
        
        filtered_taxonomic_dataframes[level] = level_filtered_taxonomic_dataframe

        # Identify which samples are below/equal to the threshold & log them
        samples_lower_equal_to_cutoff = [i for i, v in enumerate(read_counts) if v <= sample_read_counts_cutoff]
        if samples_lower_equal_to_cutoff:
            logger_output = ', '.join(samples_lower_equal_to_cutoff)
            logger.info(f"Following {level} samples had a total read count below/equal to the cutoff & were removed: {logger_output}.")
        else:
            logger.info(f"No samples were below/equal to the read counts cutoff & removed.")

    return filtered_taxonomic_dataframes


def rename_taxa(taxonomic_dfs):
    """
    Rename taxa in the taxonomic DataFrames by applying alterations, specific alterations, and homosynonyms.

    Args:
        taxonomic_dfs (dict): A dictionary with keys as taxonomic levels and values as the corresponding DataFrames.

    Returns:
        dict: A dictionary with keys as taxonomic levels and values as the renamed DataFrames.
    """
    logger.info("Renaming taxa for database compatibility.")

    resources_dir = os.path.join(os.path.dirname(__file__), 'resources')
    renaming_json_path = os.path.join(resources_dir, 'renaming.json')

    # Read the dictionaries from the JSON file
    with open(renaming_json_path, 'r') as f:
        loaded_dicts = json.load(f)

    # Access the dictionaries
    alterations, specific_alterations, homosynonyms = loaded_dicts

    renamed_dfs = {}

    for level, df in taxonomic_dfs.items():
        renamed_df = df.copy()

        # Apply alterations
        for pattern in alterations:
            renamed_df.index = renamed_df.index.str.replace(pattern, '', regex=True)

        # Apply specific alterations
        for pattern, replacement in specific_alterations.items():
            renamed_df.index = renamed_df.index.str.replace(pattern, replacement, regex=True)

        # Apply homosynonyms
        for pattern, replacement in homosynonyms.items():
            renamed_df.index = renamed_df.index.str.replace(pattern, replacement, regex=True)

        # Add the renamed DataFrame to the dictionary
        renamed_dfs[level] = renamed_df
        
        # Identify & log replaced entries
        replaced_entries = df.index[df.index != renamed_df.index]
        replacement_pairs = []
        for original_taxa_name, replacement in zip(replaced_entries, renamed_df.index[df.index != renamed_df.index]):
            replacement_pairs.append(f"Original taxa name:{original_taxa_name} - Replacement: {replacement}")
        
        if replacement_pairs:
            logger_output = ', '.join(replacement_pairs)
            logger.info(f"Original {level} taxa name(s) with their replacement(s): {logger_output}.")
        else:
            logger.info(f"No taxa namings were replaced.")

    return renamed_dfs

def check_presence_in_modelDatabase(dataframes, whichModelDatabase="full_db"):
    """
    Check if entries from the input DataFrames are in the model-database DataFrame under the same level column.
    Split the input DataFrames into two DataFrames: present and absent. 
    Add "pan" prefix to the index of the present DataFrame if the level is "Species".

    Args:
        dataframes (dict):           A dictionary containing the input DataFrames to be 
                                     checked against the model-database (currently 
                                     AGORA2 or APOLLO, or combination of both).
        whichModelDatabase (string): A string defining if AGORA2, APOLLO, or 
                                     combination of both should be used as model
                                     database to check presence in. 
                                     Allowed Input (case-insensitive): "AGORA2", "APOLLO", "full_db".
                                     Default: "full_db".

    Returns:
        dict: A dictionary containing the present and absent DataFrames for each taxonomic level.
    """
    logger.info('Checking presence of taxa in model database.')

    resources_dir = os.path.join(os.path.dirname(__file__), 'resources')
    modelDatabase_path = os.path.join(resources_dir, 'AGORA2_APOLLO_28112024.parquet')
    
    # Read in model-database as dataframe
    modelDatabase_df = pd.read_parquet(modelDatabase_path)

    # Check, if user wants to use the AGORA2 or APOLLO database, or combination of both &
    # create subset accordingly
    if whichModelDatabase.lower() == "agora2":
        updatedModelDatabase_df = modelDatabase_df[modelDatabase_df['Resource'] == 'AGORA2'].drop('Resource', axis=1)
    elif whichModelDatabase.lower() == "apollo":
        updatedModelDatabase_df = modelDatabase_df[modelDatabase_df['Resource'] == 'APOLLO'].drop('Resource', axis=1)
    else:
        updatedModelDatabase_df = modelDatabase_df.drop('Resource', axis=1)

    present_dataframes, absent_dataframes = {}, {}
    for level, input_df in dataframes.items():
        # Remove "_" from the index of the input DataFrame and find entries present in AGORA2/APOLLO
        present_mask = input_df.index.str.replace('_', ' ').isin(updatedModelDatabase_df[level])
        present_df = input_df.loc[present_mask]

        # Add "pan" prefix to the index of the present DataFrame if the level is "Species"
        if level == 'Species':
            present_df.index = 'pan' + present_df.index

        # Find entries absent in AGORA2/APOLLO
        absent_mask = ~present_mask
        absent_df = input_df.loc[absent_mask]

        present_dataframes[level] = present_df
        absent_dataframes[level] = absent_df

    return present_dataframes, absent_dataframes


def calculate_metrics(dataframes, group=None, cutoff=None, pre_mapping_read_counts=None):
    """
    Calculate & then summarize alpha diversity metrics, read counts & Firmicutes to Bacteroidetes ratio (for pyhlum)
    to compare pre-mapping & post-mapping status & evaluate the mapping coverage.

    Args:
        dataframes (dict): A dictionary with keys as taxonomic levels and values as the corresponding DataFrames.
        group (dict):      Stratified dictionary with keys as taxonomic levels and values as the corresponding DataFrames.
        cutoff (float, optional): A cut-off value for filtering out low abundance taxa. Defaults to None.
        pre_mapping_read_counts (int64 list, optional): A list containing per-sample total read counts pre-mapping,
                                allowing for taxa abundance normalization against pre-mapped total read counts.
                                Defaults to None.

    Returns:
        dicts: Dictionaries with keys as taxonomic levels and values as the calculated metrics (metrics, abundance metrics, coverage summary statistics).
    """

    metrics = {}
    abundance_metrics = {}
    summ_stats = {}

    for level, df in dataframes.items():
        if group is not None:
            df = df[group]
        
        # Group by index and sum the rows with the same name
        grouped_df = df.groupby(df.index.name).sum()

        # Get total number of taxa
        num_taxa = grouped_df.shape[0]

        # Estimate total number of named taxa, by checking for presence of "-"
        # or multiple uppercase letters in a row in taxa name
        taxa_name_contains_dash = grouped_df.index.str.contains('-')
        taxa_name_contains_uppercase = grouped_df.index.str.contains(r'[A-Z]{2,}')
        est_taxa_wo_standard_name = grouped_df[taxa_name_contains_dash | taxa_name_contains_uppercase]
        est_num_taxa_w_standard_name = len(grouped_df.index) - len(est_taxa_wo_standard_name)

        # Calculate read counts for dataframe subset & if provided, grab total read counts pre-mapping
        if pre_mapping_read_counts is not None:
            read_counts = pre_mapping_read_counts[level].sum()

            read_counts_df = grouped_df.sum()
            mean_read_counts_df = np.mean(read_counts_df)
            std_read_counts_df = np.std(read_counts_df)
        else:
            read_counts = grouped_df.sum()
            
            read_counts_df = grouped_df.sum()
            mean_read_counts_df = np.mean(read_counts)
            std_read_counts_df = np.std(read_counts)
        
        # Normalize read counts to relative abundances per taxa
        rel_abundances = grouped_df.div(read_counts)

        # Optionally apply cut-off for low abundant taxa
        if cutoff is not None:
            rel_abundances[rel_abundances <= cutoff] = 0

        # Calculate alpha diversity using the Shannon index
        shannon_index = -1 * (rel_abundances * rel_abundances.apply(np.log)).sum()
        mean_shannon_index = np.mean(shannon_index)
        std_shannon_index = np.std(shannon_index)

        # Calculate non-zero entries per column to get species richness per sample
        species_richness = (rel_abundances != 0).sum()
        mean_species_richness = np.mean(species_richness)
        std_species_richness = np.std(species_richness)

        level_metrics = {
            'read_counts': read_counts_df,
            'shannon_index': shannon_index,
        }
        
        # Calculate mean + SD, min & max relative abundance of all taxa & in 
        # how many samples a taxa is present
        level_abundance_metrics = pd.DataFrame({
            'mean': rel_abundances.mean(axis=1),
            'SD': rel_abundances.std(axis=1),
            'minimum': rel_abundances.min(axis=1),
            'maximum': rel_abundances.max(axis=1),
            'non_zero_count': (rel_abundances != 0).sum(axis=1)
        })
        
        abundance_metrics[level] = level_abundance_metrics

        if level == 'Phylum':
            # Calculate Firmicutes to Bacteroidetes ratio
            firmicutes = rel_abundances.loc['Firmicutes'] if 'Firmicutes' in rel_abundances.index else 0
            bacteroidetes = rel_abundances.loc['Bacteroidetes'] if 'Bacteroidetes' in rel_abundances.index else 0

            try:
                fb_ratio = firmicutes / bacteroidetes
            except ZeroDivisionError:
                fb_ratio = 0

            level_metrics.update({
                'firmicutes_bacteroidetes_ratio': fb_ratio,
            })

        # Add the metrics to the main dictionary
        metrics[level] = pd.DataFrame.from_dict(level_metrics)
        summ_stats[level] = [num_taxa, est_num_taxa_w_standard_name, \
                             mean_species_richness, std_species_richness, \
                             mean_shannon_index, std_shannon_index, \
                             mean_read_counts_df, std_read_counts_df]

    return metrics, abundance_metrics, summ_stats
