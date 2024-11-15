from MARS.utils import merge_files, normalize_dataframes, save_dataframes, combine_metrics
from MARS.operations import remove_clades_from_taxaNames, split_taxonomic_groups, rename_taxa, calculate_metrics, check_presence_in_agora2
import pandas as pd
import os

def process_microbial_abundances(input_file1, input_file2, output_path=None, cutoff=None, output_format="csv", stratification_file=None, flagLoneSpecies=False, taxaSplit="; ", removeCladeExtensionsFromTaxa=True):
    print(taxaSplit, flagLoneSpecies, cutoff)
    merged_dataframe = merge_files(input_file1, input_file2)
    # Added new function "remove_clades_from_taxaNames" in operations.py which can be executed optionally (default: True) - JW
    if removeCladeExtensionsFromTaxa == True:
        merged_dataframe = remove_clades_from_taxaNames(merged_dataframe)
    taxonomic_dataframes = split_taxonomic_groups(merged_dataframe, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)

    renamed_dataframes = rename_taxa(taxonomic_dataframes)
    present_dataframes, absent_dataframes = check_presence_in_agora2(renamed_dataframes)
    normalized_dataframes = normalize_dataframes(renamed_dataframes, cutoff=cutoff)
    normalized_present_dataframes, normalized_absent_dataframes = normalize_dataframes(present_dataframes, cutoff=cutoff), normalize_dataframes(absent_dataframes, cutoff=cutoff)

    pre_agora2_check_metrics = calculate_metrics(renamed_dataframes)
    post_agora2_check_metrics = calculate_metrics(present_dataframes)

    combined_metrics = combine_metrics(pre_agora2_check_metrics, post_agora2_check_metrics)

    stratification_groups = {}
    if stratification_file is not None:
        stratification = pd.read_csv(stratification_file)
        for group in stratification["group"].unique():
            # Select the columns for this group
            group_columns = list(stratification[stratification["group"] == group]["samples"])
            pre_group_metrics = calculate_metrics(renamed_dataframes, group=group_columns)
            post_group_metrics = calculate_metrics(present_dataframes, group=group_columns)

            combined_group_metrics = combine_metrics(pre_group_metrics, post_group_metrics)
            group_name = f"{group.lower()}_metrics"
            stratification_groups[group_name] = combined_group_metrics

    dataframe_groups = {'normalized': normalized_dataframes, 
                        'present': normalized_present_dataframes, 
                        'absent': normalized_absent_dataframes,
                        'metrics': combined_metrics}
    
    dataframe_groups.update(stratification_groups)

    # Save the resulting DataFrames if output_path is provided
    if output_path is not None:
        save_dataframes(dataframe_groups, output_path, output_format)

        merged_dataframe.to_csv(os.path.join(output_path,'mergedInput.csv'))

    return dataframe_groups
