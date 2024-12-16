from MARS.utils import merge_files, normalize_dataframes, normalize_dataframe, save_dataframes, combine_metrics
from MARS.operations import check_df_absolute_or_relative_counts, remove_clades_from_taxaNames, split_taxonomic_groups, rename_taxa, calculate_metrics, check_presence_in_modelDatabase, filter_samples_low_read_counts
from logging_config import setup_logger
import pandas as pd
import os
import logging


def process_microbial_abundances(input_file1, input_file2, output_path=None, cutoff=0.000001, output_format="csv", stratification_file=None, flagLoneSpecies=False, taxaSplit="; ", removeCladeExtensionsFromTaxa=True, whichModelDatabase="full_db", userDatabase_path="", sample_read_counts_cutoff=1):
    # Initialize logger to generate a MARS log file    
    logger = setup_logger('main', os.path.join(output_path, 'MARS.log'))
    logger_taxa_below_cutoff = setup_logger('taxa_below_cutoff', os.path.join(output_path, 'MARS_taxaBelowCutoff.log'))
    
    logger.info(f'INPUTS - taxaSplit: {taxaSplit}, flagLoneSpecies: {flagLoneSpecies}, sample_read_counts_cutoff: {sample_read_counts_cutoff}, cutoff: {cutoff}, removeCladeExtensionsFromTaxa: {removeCladeExtensionsFromTaxa}, whichModelDatabase: {whichModelDatabase}, userDatabase_path (if whichModelDatabase is set to "user_db"): {userDatabase_path}.')
    
    # Run MARS
    # Step 1.1: Read in input file(s) & merge dataframes when abundances + taxa names are stored seperately
    merged_dataframe = merge_files(input_file1, input_file2)
    # Step 1.2: Sum read counts (or relative abundances) for same taxa (sanity check)
    uniqueTaxa_dataframe = merged_dataframe.groupby(merged_dataframe.index.name).sum()
    # Step 1.3: Get subset of dataframe which contains species
    uniqueSpecies_dataframe = uniqueTaxa_dataframe[uniqueTaxa_dataframe.index.str.contains("s__")]
    # Step 1.4: Replace potential "_" between Genus & epithet in species name by " "
    uniqueSpecies_dataframe.index = uniqueSpecies_dataframe.index.str.replace(r'(?<!_)_(?!_)', ' ', regex=True)
    # Step 1.5: Check whether the dataframe contains absolute read counts or relative abundances and 
    # set boolean accordingly for subsequent steps
    dfvalues_are_rel_abundances = check_df_absolute_or_relative_counts(uniqueSpecies_dataframe)

    # Optional Step: Remove potential clade extensions (e.g. "clade A"; " A") from taxa namings if set true
    if removeCladeExtensionsFromTaxa == True:
        uniqueSpecies_dataframe = remove_clades_from_taxaNames(uniqueSpecies_dataframe, taxaSplit=taxaSplit)
        uniqueSpecies_dataframe = uniqueSpecies_dataframe.set_index('Taxon')

    # Step 2: Rename taxa according to resources/renaming.json to share same nomenclature as the model-database
    renamed_dataframe = rename_taxa(uniqueSpecies_dataframe)

    # Step 3: Filter out samples with too few total reads from subsequent analysis in case input 
    # is in absolute counts (not relative abundance)
    if dfvalues_are_rel_abundances == False:
        renamed_dataframe = filter_samples_low_read_counts(renamed_dataframe, sample_read_counts_cutoff=sample_read_counts_cutoff)
    
    # Step 4: Normalize & apply taxa relative abundance cutoff to renamed_dataframe
    [dataframe_afterCutoff, normalized_dataframe_afterCutoff] = normalize_dataframe(renamed_dataframe, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=cutoff)

    # Step 5: Check for presence of input taxa in a specified model database (AGORA2, APOLLO, 
    # combination of both or user-defined)
    present_dataframe, absent_dataframe = check_presence_in_modelDatabase(dataframe_afterCutoff, whichModelDatabase=whichModelDatabase, userDatabase_path=userDatabase_path, taxaSplit=taxaSplit)

    # Step 6.1: Normalize present_dataframe & absent_dataframe (with cutoff = 0,
    # because cutoff was already applied on dataframe_afterCutoff)
    [present_dataframe_afterCutoff, normalized_present_dataframe] = normalize_dataframe(present_dataframe, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=0, dataframe_to_normalize_to=renamed_dataframe)
    [absent_dataframe_afterCutoff, normalized_absent_dataframe] = normalize_dataframe(absent_dataframe, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=0, dataframe_to_normalize_to=renamed_dataframe)
    # Step 6.2: Additionally normalize present_dataframe to its own total read count to be valid input for modelling 
    # (with relative abundances of present species per sample summing to 1)
    [_, normalized_present_dataframe_adjForModelling] = normalize_dataframe(present_dataframe, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=0)

    # Step 7: Seperate merged dataframe by taxonomic levels (one df per taxonomic level)
    preMapped_dataframes = split_taxonomic_groups(dataframe_afterCutoff, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    present_dataframes = split_taxonomic_groups(present_dataframe_afterCutoff, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    absent_dataframes = split_taxonomic_groups(absent_dataframe_afterCutoff, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    preMapped_dataframes_normalized = split_taxonomic_groups(normalized_dataframe_afterCutoff, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    present_dataframes_normalized = split_taxonomic_groups(normalized_present_dataframe, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    absent_dataframes_normalized = split_taxonomic_groups(normalized_absent_dataframe, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    present_dataframes_adjForModelling = split_taxonomic_groups(normalized_present_dataframe_adjForModelling, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)

    # Step 8.1: Calculate metrics on mapping coverage & microbiome composition
    logger.info('Calculating metrices for pre-mapping dataframes.')
    pre_mapping_metrics, pre_mapping_abundance_metrics, pre_mapping_beta_diversity, pre_mapping_summ_stats = calculate_metrics(preMapped_dataframes_normalized, preMapped_dataframes)
    logger.info('Normalizing post-mapping present taxa dataframes.')
    present_post_mapping_metrics, present_post_mapping_abundance_metrics, present_post_mapping_beta_diversity, present_post_mapping_summ_stats = calculate_metrics(present_dataframes_normalized, present_dataframes)
    logger.info('Normalizing post-mapping absent taxa dataframes.')
    absent_post_mapping_metrics, absent_post_mapping_abundance_metrics, absent_post_mapping_beta_diversity, absent_post_mapping_summ_stats = calculate_metrics(absent_dataframes_normalized, absent_dataframes)
    
    # Step 8.2: Combine pre- and postMapping information of metrices, where needed
    logger.info('Combining metrics dataframes.')
    combined_metrics = combine_metrics(pre_mapping_metrics, present_post_mapping_metrics, df_type="metrics", dfvalues_are_rel_abundances=dfvalues_are_rel_abundances)
    combined_summ_stats = combine_metrics(pre_mapping_summ_stats, present_post_mapping_summ_stats, df_type="summ_stats", dfvalues_are_rel_abundances=dfvalues_are_rel_abundances)
    
    # Optional Step: If stratification groups are provided, stratify dataframe on these 
    # groups & calculate metrices on them too
    stratification_groups = {}
    stratification_groupnames = []
    if stratification_file is not None:
        stratification = pd.read_csv(stratification_file)
        for group in stratification["group"].unique():
            # Select the columns for this group
            group_columns = list(stratification[stratification["group"] == group]["samples"])
            pre_group_metrics, pre_group_abundance_metrics, pre_group_summ_stats = calculate_metrics(preMapped_dataframes_normalized, preMapped_dataframes, group=group_columns)
            post_group_metrics, post_group_abundance_metrics, post_group_summ_stats = calculate_metrics(present_dataframes_normalized, present_dataframes, group=group_columns)

            combined_group_metrics = combine_metrics(pre_group_metrics, post_group_metrics, df_type="metrics")
            combined_group_summ_stats = combine_metrics(pre_group_summ_stats, post_group_summ_stats, df_type="summ_stats")
            group_name = f"{group.lower()}_stratified_metrics"
            stratification_groupnames.append(group_name)
            stratification_groups[group_name] = [combined_group_metrics, combined_group_summ_stats, \
                                                pre_group_abundance_metrics, post_group_abundance_metrics]
        logger.info(f'Stratifying taxa dataframes using following groups: {stratification_groupnames}.')

    # Step 7: Store all results dataframe in structure & save, if output-path is provided
    dataframe_groups = {'normalized': preMapped_dataframes_normalized, 
                        'present': present_dataframes_adjForModelling, 
                        'absent': absent_dataframes_normalized,
                        'metrics': [combined_metrics, combined_summ_stats, \
                                    pre_mapping_abundance_metrics, \
                                   present_post_mapping_abundance_metrics, \
                                   absent_post_mapping_abundance_metrics, \
                                   pre_mapping_beta_diversity, \
                                   present_post_mapping_beta_diversity, \
                                   absent_post_mapping_beta_diversity]}
    
    dataframe_groups.update(stratification_groups)

    if output_path is not None:
        logger.info(f'Saving output to {output_path}.')
        save_dataframes(dataframe_groups, output_path, output_format)

        merged_dataframe.to_csv(os.path.join(output_path,'mergedInput.csv'))

    return dataframe_groups
