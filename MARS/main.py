from MARS.utils import merge_files, normalize_dataframes, save_dataframes, combine_metrics
from MARS.operations import check_df_absolute_or_relative_counts, remove_clades_from_taxaNames, split_taxonomic_groups, rename_taxa, calculate_metrics, check_presence_in_modelDatabase, filter_samples_low_read_counts
from logging_config import setup_logger
import pandas as pd
import os
import logging


def process_microbial_abundances(input_file1, input_file2, output_path=None, cutoff=None, output_format="csv", stratification_file=None, flagLoneSpecies=False, taxaSplit="; ", removeCladeExtensionsFromTaxa=True, whichModelDatabase="full_db", userDatabase_path=""):
    # Initialize logger to generate a MARS log file    
    logger = setup_logger('main', os.path.join(output_path, 'MARS.log'))
    logger_taxa_below_cutoff = setup_logger('taxa_below_cutoff', os.path.join(output_path, 'MARS_taxaBelowCutoff.log'))
    
    logger.info(f'INPUTS - taxaSplit: {taxaSplit}, flagLoneSpecies: {flagLoneSpecies}, cutoff: {cutoff}, removeCladeExtensionsFromTaxa: {removeCladeExtensionsFromTaxa}, whichModelDatabase: {whichModelDatabase}.')
    
    # Run MARS
    # Step 1.1: Read in input file(s) & merge dataframes when abundances + taxa names are stored seperately
    merged_dataframe = merge_files(input_file1, input_file2)
    # Step 1.2: Replace pot. "_" between Genus & epithet in species name by " "
    merged_dataframe.index = merged_dataframe.index.str.replace(r'(?<!_)_(?!_)', ' ', regex=True)
    # Step 1.3: Check whether dataframe contains absolute read counts or relative abundances and set boolean for acording subsequent steps
    dfvalues_are_rel_abundances = check_df_absolute_or_relative_counts(merged_dataframe)

    # Optional Step: Remove pot. clade extensions (e.g. clade A; A) from taxa namings if set true
    if removeCladeExtensionsFromTaxa == True:
        merged_dataframe = remove_clades_from_taxaNames(merged_dataframe, taxaSplit=taxaSplit)
    
    # Step 2: Seperate merged dataframe by taxonomic levels (one df per taxonomic level)
    taxonomic_dataframes = split_taxonomic_groups(merged_dataframe, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    
    # Optional Step: Filter out samples with too view total reads from subsequent analysis in case input is in absolute counts (not relative abundance)
    if dfvalues_are_rel_abundances == False:
        taxonomic_dataframes = filter_samples_low_read_counts(taxonomic_dataframes, sample_read_counts_cutoff=1)
    
    # Step 3: Rename taxa according to resources/renaming.json to share same nomenclature as the model-database
    renamed_dataframes = rename_taxa(taxonomic_dataframes)
    
    # Step 4: Check for presence of input taxa in a specified model database (AGORA2, APOLLO, combination of both or user-defined)
    present_dataframes, absent_dataframes = check_presence_in_modelDatabase(renamed_dataframes, whichModelDatabase=whichModelDatabase, userDatabase_path=userDatabase_path)
    
    # Step 5: Normalize read counts to obtain realtive abundances
    logger.info('Normalizing pre-mapping dataframes.')
    normalized_dataframes = normalize_dataframes(renamed_dataframes, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=cutoff)
    logger.info('Normalizing post-mapping present & absent taxa dataframes.')
    normalized_present_dataframes, normalized_absent_dataframes = normalize_dataframes(present_dataframes, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=cutoff), normalize_dataframes(absent_dataframes, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=cutoff, pre_mapping_read_counts=renamed_dataframes)
    
    # Step 6.1: Calculate metrics on mapping coverage & microbiome composition
    logger.info('Calculating metrices for pre-mapping dataframes.')
    pre_mapping_metrics, pre_mapping_abundance_metrics, pre_mapping_summ_stats = calculate_metrics(renamed_dataframes, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=cutoff)
    logger.info('Normalizing post-mapping present taxa dataframes.')
    present_post_mapping_metrics, present_post_mapping_abundance_metrics, present_post_mapping_summ_stats = calculate_metrics(present_dataframes, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=cutoff, pre_mapping_read_counts=renamed_dataframes)
    logger.info('Normalizing post-mapping absent taxa dataframes.')
    absent_post_mapping_metrics, absent_post_mapping_abundance_metrics, absent_post_mapping_summ_stats = calculate_metrics(absent_dataframes, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=cutoff, pre_mapping_read_counts=renamed_dataframes)
    
    # Step 6.2: Combine pre- and postMapping information of metrices, where needed
    logger.info('Combining metrics dataframes.')
    combined_metrics = combine_metrics(pre_mapping_metrics, present_post_mapping_metrics, df_type="metrics", dfvalues_are_rel_abundances=dfvalues_are_rel_abundances)
    combined_summ_stats = combine_metrics(pre_mapping_summ_stats, present_post_mapping_summ_stats, df_type="summ_stats", dfvalues_are_rel_abundances=dfvalues_are_rel_abundances)
    
    # Optional Step: If stratification groups are provided, stratify dataframe on these groups & calculate metrices on them too
    stratification_groups = {}
    stratification_groupnames = []
    if stratification_file is not None:
        stratification = pd.read_csv(stratification_file)
        for group in stratification["group"].unique():
            # Select the columns for this group
            group_columns = list(stratification[stratification["group"] == group]["samples"])
            pre_group_metrics, pre_group_abundance_metrics, pre_group_summ_stats = calculate_metrics(renamed_dataframes, group=group_columns, cutoff=cutoff)
            post_group_metrics, post_group_abundance_metrics, post_group_summ_stats = calculate_metrics(present_dataframes, group=group_columns, cutoff=cutoff, pre_mapping_read_counts=renamed_dataframes)

            combined_group_metrics = combine_metrics(pre_group_metrics, post_group_metrics, df_type="metrics")
            combined_group_summ_stats = combine_metrics(pre_group_summ_stats, post_group_summ_stats, df_type="summ_stats")
            group_name = f"{group.lower()}_stratified_metrics"
            stratification_groupnames.append(group_name)
            stratification_groups[group_name] = [combined_group_metrics, combined_group_summ_stats, \
                                                pre_group_abundance_metrics, post_group_abundance_metrics]
        logger.info(f'Stratifying taxa dataframes using following groups: {stratification_groupnames}.')
    
    # Step 7: Store all results dataframe in structure & save, if output-path is provided
    dataframe_groups = {'normalized': normalized_dataframes, 
                        'present': normalized_present_dataframes, 
                        'absent': normalized_absent_dataframes,
                        'metrics': [combined_metrics, combined_summ_stats, \
                                    pre_mapping_abundance_metrics, \
                                   present_post_mapping_abundance_metrics, \
                                   absent_post_mapping_abundance_metrics]}
    
    dataframe_groups.update(stratification_groups)

    if output_path is not None:
        logger.info(f'Saving output to {output_path}.')
        save_dataframes(dataframe_groups, output_path, output_format)

        merged_dataframe.to_csv(os.path.join(output_path,'mergedInput.csv'))

    return dataframe_groups
