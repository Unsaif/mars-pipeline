from MARS.utils import merge_files, normalize_dataframes, save_dataframes, combine_metrics
from MARS.operations import remove_clades_from_taxaNames, split_taxonomic_groups, rename_taxa, calculate_metrics, check_presence_in_modelDatabase
from logging_config import setup_logger
import pandas as pd
import os
import logging


def process_microbial_abundances(input_file1, input_file2, output_path=None, cutoff=None, output_format="csv", stratification_file=None, flagLoneSpecies=False, taxaSplit="; ", removeCladeExtensionsFromTaxa=True, whichModelDatabase="both"):
    # Initialize logger to generate a MARS log file
    logging_file_name = 'MARS.log'
    logging_file_path = os.path.join(output_path, logging_file_name)
    
    logger = setup_logger('main', logging_file_path)
    
    logger.info(f'taxaSplit: {taxaSplit}, flagLoneSpecies: {flagLoneSpecies}, cutoff: {cutoff}, removeCladeExtensionsFromTaxa: {removeCladeExtensionsFromTaxa}, whichModelDatabase: {whichModelDatabase}.')
    
    # Run MARS
    merged_dataframe = merge_files(input_file1, input_file2)

    if removeCladeExtensionsFromTaxa == True:
        merged_dataframe = remove_clades_from_taxaNames(merged_dataframe, taxaSplit=taxaSplit)

    taxonomic_dataframes = split_taxonomic_groups(merged_dataframe, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    # filtered_taxonomic_dataframes = filter_samples_low_read_counts(taxonomic_dataframes)

    renamed_dataframes = rename_taxa(taxonomic_dataframes)

    present_dataframes, absent_dataframes = check_presence_in_modelDatabase(renamed_dataframes, whichModelDatabase=whichModelDatabase)
    
    logger.info('Normalizing pre-mapping dataframes.')
    normalized_dataframes = normalize_dataframes(renamed_dataframes, cutoff=cutoff)
    logger.info('Normalizing post-mapping present & absent taxa dataframes.')
    normalized_present_dataframes, normalized_absent_dataframes = normalize_dataframes(present_dataframes, cutoff=cutoff), normalize_dataframes(absent_dataframes, cutoff=cutoff, pre_mapping_read_counts=renamed_dataframes)
    
    logger.info('Calculating metrices for pre-mapping dataframes.')
    pre_agora2_check_metrics, pre_mapping_abundance_metrics, pre_mapping_summ_stats = calculate_metrics(renamed_dataframes, cutoff=cutoff)
    logger.info('Normalizing post-mapping present taxa dataframes.')
    post_agora2_check_metrics, present_post_mapping_abundance_metrics, present_post_mapping_summ_stats = calculate_metrics(present_dataframes, cutoff=cutoff, pre_mapping_read_counts=renamed_dataframes)
    logger.info('Normalizing post-mapping absent taxa dataframes.')
    absent_post_mapping_metrics, absent_post_mapping_abundance_metrics, absent_post_mapping_summ_stats = calculate_metrics(absent_dataframes, cutoff=cutoff, pre_mapping_read_counts=renamed_dataframes)
    
    logger.info('Combining metrics dataframes.')
    combined_metrics = combine_metrics(pre_agora2_check_metrics, post_agora2_check_metrics, df_type="metrics")
    combined_summ_stats = combine_metrics(pre_mapping_summ_stats, present_post_mapping_summ_stats, df_type="summ_stats")

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

    dataframe_groups = {'normalized': normalized_dataframes, 
                        'present': normalized_present_dataframes, 
                        'absent': normalized_absent_dataframes,
                        'metrics': [combined_metrics, combined_summ_stats, \
                                    pre_mapping_abundance_metrics, \
                                   present_post_mapping_abundance_metrics, \
                                   absent_post_mapping_abundance_metrics]}
    
    dataframe_groups.update(stratification_groups)

    # Save the resulting DataFrames if output_path is provided
    if output_path is not None:
        save_dataframes(dataframe_groups, output_path, output_format)

        merged_dataframe.to_csv(os.path.join(output_path,'mergedInput.csv'))

    return dataframe_groups
