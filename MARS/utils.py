import numpy as np
import pandas as pd
import os
import logging

logger = logging.getLogger('main.utils')
logger_taxa_below_cutoff = logging.getLogger('taxa_below_cutoff.utils')

def read_file_as_dataframe(file_path, header):
    # Extract the file extension
    try:
        file_extension = file_path.type
    except AttributeError:
        ind = file_path.find('.')
        extension = file_path[ind+1:]
        if extension == 'txt':
            file_extension = 'text/plain'
        elif extension == 'csv':
            file_extension = 'text/csv'
        elif extension == 'xlsx':
            file_extension = 'spreadsheet'
        else:
            file_extension = 'not found'
    # Read the file based on its extension
    if file_extension == 'text/plain' or file_extension == 'text/tab-separated-values':
        # Assuming the txt file is delimited (e.g., CSV)
        return pd.read_csv(file_path, sep='\t', index_col=[0], low_memory=False, header=header)  # Update delimiter if necessary
    elif file_extension == 'text/csv' or file_extension == 'application/vnd.ms-excel':
        return pd.read_csv(file_path, index_col=[0], low_memory=False, header=header)
    elif file_extension == 'spreadsheet' or file_extension == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet':
        # Removed low_memory input argument as it seems it is not longer accepted in the pd.read_excel function
        # return pd.read_excel(file_path, engine='openpyxl', index_col=[0], low_memory=False, header=header)
        return pd.read_excel(file_path, engine='openpyxl', index_col=[0], header=header)
    else:
        raise ValueError(f"Unsupported file type: {file_extension}")

def merge_files(file_path1, file_path2):
    """
    Merge two files using pandas and return the merged DataFrame.

    Args:
        file_path1 (str): The file path of the first input file.
        file_path2 (str): The file path of the second input file.

    Returns:
        pd.DataFrame: The merged DataFrame.
    """

    # Read input files into pandas DataFrames  
    
    if file_path2 == None:
        df1 = read_file_as_dataframe(file_path1, 0)
        merged_df = df1
    else:
        df1 = read_file_as_dataframe(file_path1, 0)
        df2 = read_file_as_dataframe(file_path2, 1) 
        # Merge DataFrames using their index values
        merged_df = pd.merge(df1, df2, left_index=True, right_index=True, how='inner')
        # Drop the 'Confidence' column
        try:
            merged_df = merged_df.drop(columns=['Confidence'])
        except:
            pass

        # Reset the index and set the 'Taxon' column as the new index
        merged_df = merged_df.reset_index(drop=True).set_index('Taxon')

    return merged_df

def normalize_dataframes(dataframes, cutoff=None, pre_mapping_read_counts=None):
    """
    Normalize the taxonomic DataFrames by grouping and summing rows with the same name,
    and then normalizing each column so that the sum of each column is 1. Optionally, a
    cut-off can be provided to filter out low abundance taxa before normalization.

    Args:
        dataframes (dict): A dictionary with keys as taxonomic levels and values as the corresponding DataFrames.
        cutoff (float, optional): A cut-off value for filtering out low abundance taxa. Defaults to None.
        pre_mapping_read_counts (int64 list, optional): A list containing per-sample total read counts pre-mapping,
                                allowing for taxa abundance normalization against pre-mapped total read counts.
                                Defaults to None.

    Returns:
        dict: A dictionary with keys as taxonomic levels and values as the normalized DataFrames.
    """

    normalized_dfs = {}

    for level, df in dataframes.items():
        # Group by index and sum the rows with the same name
        grouped_df = df.groupby(df.index.name).sum()

        # Normalize each column so that the sum of each column is 1 (either
        # to pre-mapped total read counts, or to the subset read counts)
        if pre_mapping_read_counts is not None:
            read_counts = pre_mapping_read_counts[level].sum()
        else:
            read_counts = grouped_df.sum()
        
        # Normalize read counts to get relative abundances of taxa
        rel_abundances_df = grouped_df.div(read_counts)
  
        # Optionally apply cut-off for low abundance taxa. Coincidentally
        # also fixes empty cells to 0s.
        if cutoff is not None:
            rel_abundances_df[rel_abundances_df <= cutoff] = 0

            # Identify which taxa in which samples are below cutoff threshold & set to 0, log them
            entries_below_cutoff = rel_abundances_df[rel_abundances_df <= cutoff].stack().index.tolist()

            if entries_below_cutoff:
                logger.info(f"{level} taxa were below the cutoff & are listed in seperate log-file.")
                logger_taxa_below_cutoff.info(f"{level} taxa whose rel.abundance was below the cutoff & therefore set to 0: {entries_below_cutoff}")
            else:
                logger.info(f"No {level} taxa were below the cutoff.")
            
        # Add the normalized DataFrame to the dictionary
        normalized_dfs[level] = rel_abundances_df
        
    return normalized_dfs


def combine_metrics(metrics1, metrics2, df_type="metrics"):
    """
    Combine the metrics from two different sets of taxonomic DataFrames into a single DataFrame for each level.

    Args:
        metrics1 (dict): A dictionary with keys as taxonomic levels and values as the calculated metrics for the first group.
        metrics2 (dict): A dictionary with keys as taxonomic levels and values as the calculated metrics for the second group.

    Returns:
        dict: A dictionary with keys as taxonomic levels and values as the combined DataFrames.
    """
    combined_metrics = {}

    for level in metrics1.keys():
        level_metrics_pre_mapping = metrics1[level]
        level_metrics_post_mapping = metrics2[level]

        if df_type == "metrics":
            combined_level_metrics = {}

            for metric_name in level_metrics_pre_mapping.keys():
                
                combined_metric = pd.DataFrame([level_metrics_pre_mapping[metric_name], level_metrics_post_mapping[metric_name]])
                combined_metric.index = ['pre AGORA2 mapping', 'post AGORA2 mapping']
                combined_level_metrics[metric_name] = combined_metric

            combined_metrics[level] = combined_level_metrics

        elif df_type == "summ_stats":
            combined_metric = pd.DataFrame({'Pre mapping': level_metrics_pre_mapping, \
                                            'Post mapping': level_metrics_post_mapping})
            
            combined_metric['Mapping coverage'] = np.nan
            combined_metric.loc[0, 'Mapping coverage'] = combined_metric.loc[0, 'Post mapping'] / combined_metric.loc[0, 'Pre mapping']
            combined_metric.loc[2, 'Mapping coverage'] = combined_metric.loc[2, 'Post mapping'] / combined_metric.loc[2, 'Pre mapping']
            combined_metric.loc[4, 'Mapping coverage'] = combined_metric.loc[4, 'Post mapping'] - combined_metric.loc[4, 'Pre mapping']
            combined_metric.loc[6, 'Mapping coverage'] = combined_metric.loc[6, 'Post mapping'] / combined_metric.loc[6, 'Pre mapping']
            
            combined_metric['Description'] = ['The number of taxa across all samples. MappingCoverage = post mapping/pre mapping.', \
                                              'The estimated number of taxa following standard nomenclature. Estimated by excluding all taxa whose names contain "-" &/or multiple uppercase letters in a row.', \
                                              'Mean number of species across samples. MappingCoverage = post mapping/pre mapping.', \
                                              'Standard deviation of species richness.', 'Mean alpha-diversity (calc. by shannon index) across samples. MappingCoverage = post mapping - pre mapping.', \
                                             'Standard deviation of shannon index.', 'Mean number of reads across samples. MappingCoverage = post mapping/pre mapping.', \
                                             'Standard deviation of read counts.']

            combined_metric.index = ['Total number of taxa', 'Estimated total number of named taxa', 'Mean species richness', \
                                     'Std species richness', 'Mean shannon index', 'Std shannon index', \
                                     'Mean read counts', 'Std read counts']
            
            combined_metrics[level] = combined_metric

    return combined_metrics


def save_dataframes(dataframe_groups, output_path, output_format):
    """
    Save DataFrames to the specified output path in the given format.

    Args:
        dataframe_groups (dict): A dictionary of dictionaries containing DataFrames.
                                  Outer dictionary keys represent group names, and inner dictionaries contain
                                  DataFrames for each taxonomic level.
        output_path (str): The directory path where the output files will be saved.
        output_format (str): The format to save the DataFrames ('csv', 'excel', 'parquet', or 'json').
    """
    os.makedirs(output_path, exist_ok=True)

    for group_name, values in dataframe_groups.items():
        group_output_path = os.path.join(output_path, group_name)
        os.makedirs(group_output_path, exist_ok=True)
        if group_name == "metrics":
            for level, metrics_dataframes in values[0].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                for metric_name, df in metrics_dataframes.items():
                    file_name = f"{metric_name}.{output_format}"
                    file_path = os.path.join(level_output_path, file_name)

                    if output_format == "csv":
                        df.to_csv(file_path)
                    elif output_format == "txt" or output_format == "tsv":
                        df.to_csv(file_path, sep='\t')
                    elif output_format == "excel":
                        df.to_excel(file_path)
                    elif output_format == "parquet":
                        df.to_parquet(file_path)
                    elif output_format == "json":
                        df.to_json(file_path)
                    else:
                        raise ValueError(f"Unsupported output format: {output_format}")

            for level, summ_stats_df in values[1].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"summ_stats_{level}.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    summ_stats_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    summ_stats_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    summ_stats_df.to_excel(file_path)
                elif output_format == "parquet":
                    summ_stats_df.to_parquet(file_path)
                elif output_format == "json":
                    summ_stats_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")

            for level, abundance_metrics_df in values[2].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"preMapping_abundanceMetrics_{level}.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    abundance_metrics_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    abundance_metrics_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    abundance_metrics_df.to_excel(file_path)
                elif output_format == "parquet":
                    abundance_metrics_df.to_parquet(file_path)
                elif output_format == "json":
                    abundance_metrics_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")

            for level, abundance_metrics_df in values[3].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"postMapping_present_abundanceMetrics_{level}.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    abundance_metrics_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    abundance_metrics_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    abundance_metrics_df.to_excel(file_path)
                elif output_format == "parquet":
                    abundance_metrics_df.to_parquet(file_path)
                elif output_format == "json":
                    abundance_metrics_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")
    
            for level, abundance_metrics_df in values[4].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"postMapping_absent_abundanceMetrics_{level}.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    abundance_metrics_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    abundance_metrics_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    abundance_metrics_df.to_excel(file_path)
                elif output_format == "parquet":
                    abundance_metrics_df.to_parquet(file_path)
                elif output_format == "json":
                    abundance_metrics_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")
        
        elif "_stratified_metrics" in group_name:
            for level, metrics_dataframes in values[0].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                for metric_name, df in metrics_dataframes.items():
                    file_name = f"{metric_name}_stratified.{output_format}"
                    file_path = os.path.join(level_output_path, file_name)

                    if output_format == "csv":
                        df.to_csv(file_path)
                    elif output_format == "txt" or output_format == "tsv":
                        df.to_csv(file_path, sep='\t')
                    elif output_format == "excel":
                        df.to_excel(file_path)
                    elif output_format == "parquet":
                        df.to_parquet(file_path)
                    elif output_format == "json":
                        df.to_json(file_path)
                    else:
                        raise ValueError(f"Unsupported output format: {output_format}")

            for level, summ_stats_df in values[1].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"summ_stats_{level}_stratified.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    summ_stats_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    summ_stats_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    summ_stats_df.to_excel(file_path)
                elif output_format == "parquet":
                    summ_stats_df.to_parquet(file_path)
                elif output_format == "json":
                    summ_stats_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")

            for level, abundance_metrics_df in values[2].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"preMapping_abundanceMetrics_{level}_stratified.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    abundance_metrics_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    abundance_metrics_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    abundance_metrics_df.to_excel(file_path)
                elif output_format == "parquet":
                    abundance_metrics_df.to_parquet(file_path)
                elif output_format == "json":
                    abundance_metrics_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")

            for level, abundance_metrics_df in values[3].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"postMapping_abundanceMetrics_{level}_stratified.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    abundance_metrics_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    abundance_metrics_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    abundance_metrics_df.to_excel(file_path)
                elif output_format == "parquet":
                    abundance_metrics_df.to_parquet(file_path)
                elif output_format == "json":
                    abundance_metrics_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")
    
        else:
            for level, df in values.items():
                file_name = f"{group_name}_{level.lower()}.{output_format}"
                file_path = os.path.join(group_output_path, file_name)

                if output_format == "csv":
                    df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    df.to_excel(file_path)
                elif output_format == "parquet":
                    df.to_parquet(file_path)
                elif output_format == "json":
                    df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")
            
