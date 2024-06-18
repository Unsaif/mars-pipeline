import pandas as pd
import os

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
            file_extension == 'spreadsheet'

    # Read the file based on its extension
    if file_extension == 'text/plain' or file_extension == 'text/tab-separated-values':
        # Assuming the txt file is delimited (e.g., CSV)
        return pd.read_csv(file_path, sep='\t', index_col=[0], low_memory=False, header=header)  # Update delimiter if necessary
    elif file_extension == 'text/csv':
        return pd.read_csv(file_path, index_col=[0], low_memory=False, header=header)
    elif file_extension == 'spreadsheet':
        return pd.read_excel(file_path, engine='openpyxl', index_col=[0], low_memory=False, header=header)
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

def normalize_dataframes(dataframes, cutoff=None):
    """
    Normalize the taxonomic DataFrames by grouping and summing rows with the same name,
    and then normalizing each column so that the sum of each column is 1. Optionally, a
    cut-off can be provided to filter out low abundance taxa before normalization.

    Args:
        dataframes (dict): A dictionary with keys as taxonomic levels and values as the corresponding DataFrames.
        cutoff (float, optional): A cut-off value for filtering out low abundance taxa. Defaults to None.

    Returns:
        dict: A dictionary with keys as taxonomic levels and values as the normalized DataFrames.
    """

    normalized_dfs = {}

    for level, df in dataframes.items():
        # Group by index and sum the rows with the same name
        grouped_df = df.groupby(df.index.name).sum()

        # Optionally apply cut-off for low abundance taxa
        if cutoff is not None:
            grouped_df = grouped_df[grouped_df >= cutoff]

        # Normalize each column so that the sum of each column is 1
        normalized_df = grouped_df / grouped_df.sum()

        # Add the normalized DataFrame to the dictionary
        normalized_dfs[level] = normalized_df

    return normalized_dfs

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

    for group_name, dataframes in dataframe_groups.items():
        group_output_path = os.path.join(output_path, group_name)
        os.makedirs(group_output_path, exist_ok=True)
        if "metrics" in group_name:
            for level, metrics_dataframes in dataframes.items():
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
        else:
            for level, df in dataframes.items():
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
            
def combine_metrics(metrics1, metrics2):
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
        level_metrics1 = metrics1[level]
        level_metrics2 = metrics2[level]

        combined_level_metrics = {}

        for metric_name in level_metrics1.keys():
            combined_metric = pd.DataFrame([level_metrics1[metric_name], level_metrics2[metric_name]])
            combined_metric.index = ['pre AGORA2 mapping', 'post AGORA2 mapping']
            combined_level_metrics[metric_name] = combined_metric

        combined_metrics[level] = combined_level_metrics

    return combined_metrics