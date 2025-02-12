# app.py
import streamlit as st
import pandas as pd
import os
from io import BytesIO
from MARS.utils import normalize_dataframe, save_dataframes, combine_metrics, read_file_as_dataframe, merge_files
from MARS.operations import load_input_and_preprocess, check_df_absolute_or_relative_counts, remove_clades_from_taxaNames, concatenate_genus_and_species_names, split_taxonomic_groups, rename_taxa, calculate_metrics, check_presence_in_modelDatabase, filter_samples_low_read_counts
from ANT.ant import ant

st.set_page_config(layout="wide")

# --- Helper functions ---
def convert_df(df, file_format, index=False):
    """Converts a Pandas DataFrame to a specified file format for download."""
    if file_format == "xlsx":
        output = BytesIO()
        with pd.ExcelWriter(output, engine="openpyxl") as writer:
            df.to_excel(writer, index=index)
        data = output.getvalue()
        return data
    elif file_format == "txt":
        return df.to_csv(index=index, sep="\t").encode("utf-8")
    else:
        return df.to_csv(index=index).encode("utf-8")

def file_to_list(uploaded_file):
    """Reads a file and returns the data as a list.  Handles different file types."""
    if uploaded_file is not None:
        file_details = {
            "FileName": uploaded_file.name,
            "FileType": uploaded_file.type,
            "FileSize": uploaded_file.size,
        }
        st.write(file_details)

        if uploaded_file.type == "text/plain" or uploaded_file.type == "text/tab-separated-values":
            # For txt files, read each line into a list
            bytes_data = uploaded_file.getvalue()
            str_data = bytes_data.decode("utf-8")
            lines = str_data.split("\n")
            return [line for line in lines if line]
        elif uploaded_file.type == "text/csv":
            try:
                df = pd.read_csv(uploaded_file)
                # Assuming the species names are in the first column
                species_list = df.iloc[:, 0].tolist()
                return species_list
            except Exception as e:
                st.write("Could not read file: ", e)
                return None
        elif "spreadsheet" in uploaded_file.type:
            try:
                df = pd.read_excel(uploaded_file)
                # Assuming the species names are in the first column
                species_list = df.iloc[:, 0].tolist()
                return species_list
            except Exception as e:
                st.write("Could not read file: ", e)
                return None
        else:
            st.write("Unsupported file type: ", uploaded_file.type)
            return None
    else:
        return None

# --- Streamlit app ---
st.image("images/logo2.png", use_column_width=True)

st.title("MARS")
st.write("""
""")

# --- Input files ---
col1, col2, col3, col4 = st.columns(4)
input_file1 = col1.file_uploader("Upload Taxonomy Table", type=['csv', 'tsv', 'txt', 'xlsx'])
input_file2 = col2.file_uploader("Upload Feature Table", type=['csv', 'tsv', 'txt', 'xlsx'])
input_file3 = col3.file_uploader("Upload Combined Table", type=['csv', 'tsv', 'txt', 'xlsx'])
uploaded_resource_file = col4.file_uploader(
    "(Optional) Upload your own resource file for ANT",
    type=["txt", "csv", "xlsx"],
)

# --- Sidebar for additional parameters ---
st.sidebar.header("Parameters")
output_path = st.sidebar.text_input("Output Path (optional)", None)
cutoff = st.sidebar.slider("Cutoff", 0.0, 100.0, 0.0)
output_format = st.sidebar.radio("Output Format", ["csv", "txt", "xlsx"])
stratification_file = st.sidebar.file_uploader("Stratification File", type=['csv', 'txt', 'xlsx'])
skip_ant = st.sidebar.checkbox('Skip ANT')
flagLoneSpecies = st.sidebar.checkbox('Genus name not present in the s__ taxonomic identifier?')
taxaSplit = st.sidebar.text_input('Delimiter used to seperate taxonomic levels', ';')
removeCladeExtensionsFromTaxa = st.sidebar.checkbox('Remove Clade Extensions from Taxa Names', True)
sample_read_counts_cutoff = st.sidebar.number_input("Sample Read Counts Cutoff", min_value=1, value=1)
whichModelDatabase = st.sidebar.selectbox("Model Database", ["full_db", "AGORA2", "APOLLO", "user_db"], index=0)
userDatabase_path = st.sidebar.text_input("User Database Path (if user_db is selected)", "")


# --- Main Processing ---
if (input_file1 and input_file2) or input_file3:
    if st.button("Run MARS"):
        try:
            with st.spinner("Running MARS..."):
                # Step 1: Check input data & preprocess
                if input_file3:
                    preprocessed_dataframe, dfvalues_are_rel_abundances = load_input_and_preprocess(input_file3, taxaSplit=taxaSplit)
                    st.success("Input file loaded and preprocessed.")
                else:
                    preprocessed_dataframe, dfvalues_are_rel_abundances = load_input_and_preprocess(input_file1, input_file2, taxaSplit)
                    st.success("Input files loaded and preprocessed.")
                    
                # Optional Step: Remove potential clade extensions (e.g. "clade A"; " A") from taxa namings if set true
                if removeCladeExtensionsFromTaxa:
                    preprocessed_dataframe = remove_clades_from_taxaNames(preprocessed_dataframe, taxaSplit=taxaSplit)
                    st.info("Clade extensions removed.")

                # Optional Step: Concatenate genus name with species epithet if both are present, otherwise leave species column unchanged
                if flagLoneSpecies:
                    preprocessed_dataframe = concatenate_genus_and_species_names(preprocessed_dataframe, taxaSplit=taxaSplit)
                    st.info("Genus and species names concatenated.")

                # Step 2: Rename taxa according to resources/renaming.json to share same nomenclature as the model-data
                renamed_dataframe = rename_taxa(preprocessed_dataframe)
                st.info("Taxa renamed.")

                # Step 3: Filter out samples with too few total reads from subsequent analysis in case input is in absolute read counts (not relative abundance)
                if not dfvalues_are_rel_abundances:
                    renamed_dataframe = filter_samples_low_read_counts(renamed_dataframe, sample_read_counts_cutoff=sample_read_counts_cutoff)
                    st.info(f"Samples with read counts below {sample_read_counts_cutoff} filtered out.")


                # --- ANT Integration ---
                def replace_homosynonyms(index):
                    """
                    Replace all occurrences of old substring with new substring in dataframe index.
                    """
                    for old, new in homosynonyms.items():
                        index = index.replace(old, new)
                    return index

                if not skip_ant:
                    resource = file_to_list(uploaded_resource_file)
                    parts = renamed_dataframe.index.split(taxaSplit)
                    species = [part for part in parts if "s__" in part]
                    species = [sub.replace('_', ' ') for sub in species]

                    species_in_resource, homosynonyms, ncbi_tax_id, all_found_homosynoyms = ant(
                        species, resource
                    )
                    # Replace species names by their newly found homosynonyms 
                    if homosynonyms:
                        with st.spinner("Renaming found homosynonyms..."):
                            homosynonyms = {key.replace(' ', '_'): value.replace(' ', '_') for key, value in homosynonyms.items()}
                            renamed_dataframe.index = renamed_dataframe.index.map(replace_homosynonyms)
                    st.success("Renaming with found homosynonyms: Success!")
                    
                    # Display ANT results
                    col1, col2 = st.columns(2)
                    col1.metric("Number of species found in resource", len(species_in_resource))
                    col2.metric("Number of homosynonyms found", len(homosynonyms))

                    st.success("ANT analysis complete!")
                else:
                    st.info("ANT analysis skipped.")


                # Step 4: Normalize & apply taxa relative abundance cutoff to renamed_dataframe
                dataframe_afterCutoff, normalized_dataframe_afterCutoff = normalize_dataframe(renamed_dataframe, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=cutoff)
                st.info("Data normalized and cutoff applied.")

                # Step 5: Check for presence of input taxa in a specified model database (AGORA2, APOLLO, combination of both or user-defined)
                present_dataframe, absent_dataframe = check_presence_in_modelDatabase(dataframe_afterCutoff, whichModelDatabase=whichModelDatabase, userDatabase_path=userDatabase_path, taxaSplit=taxaSplit)
                if present_dataframe.empty:
                    st.error("No species found in the specified model database. Please check your input and database selection.")
                    raise ValueError("No species from the input data were found & could be mapped to the reconstruction database.")
                st.info("Presence in model database checked.")

                # Step 6.1: Normalize present_dataframe & absent_dataframe (with cutoff = 0, because cutoff was already applied on dataframe_afterCutoff)
                present_dataframe_afterCutoff, normalized_present_dataframe = normalize_dataframe(present_dataframe, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=0, dataframe_to_normalize_to=dataframe_afterCutoff)
                absent_dataframe_afterCutoff, normalized_absent_dataframe = normalize_dataframe(absent_dataframe, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=0, dataframe_to_normalize_to=dataframe_afterCutoff)
                # Step 6.2: Additionally normalize present_dataframe to its own total read count to be valid input for modelling (with relative abundances of present species per sample summing to 1)
                _, normalized_present_dataframe_adjForModelling = normalize_dataframe(present_dataframe, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=0)
                st.info("Dataframes normalized.")

                # Step 7.1: Seperate normalized dataframes by taxonomic levels (one df per taxonomic level, for both normalized & not-normalized dataframes)
                preMapped_dataframes = split_taxonomic_groups(dataframe_afterCutoff, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
                present_dataframes = split_taxonomic_groups(present_dataframe_afterCutoff, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
                absent_dataframes = split_taxonomic_groups(absent_dataframe_afterCutoff, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
                preMapped_dataframes_normalized = split_taxonomic_groups(normalized_dataframe_afterCutoff, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
                present_dataframes_normalized = split_taxonomic_groups(normalized_present_dataframe, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
                absent_dataframes_normalized = split_taxonomic_groups(normalized_absent_dataframe, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
                present_dataframes_adjForModelling = split_taxonomic_groups(normalized_present_dataframe_adjForModelling, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
                st.info("Taxonomic groups split.")

                # Step 7.2: Add "pan" prefix to species names (in index) of present_dataframes_adjForModelling (because pan-species reconstructions will be used in MgPipe)
                present_dataframes_adjForModelling['Species'].index = 'pan' + present_dataframes_adjForModelling['Species'].index
                st.info("\"pan\" prefix added to species names.")

                # Step 8.1: Calculate metrics on mapping coverage & microbiome composition)
                pre_mapping_metrics, pre_mapping_abundance_metrics, pre_mapping_beta_diversity, pre_mapping_summ_stats = calculate_metrics(preMapped_dataframes_normalized, preMapped_dataframes)
                present_post_mapping_metrics, present_post_mapping_abundance_metrics, present_post_mapping_beta_diversity, present_post_mapping_summ_stats = calculate_metrics(present_dataframes_normalized, present_dataframes)
                absent_post_mapping_metrics, absent_post_mapping_abundance_metrics, absent_post_mapping_beta_diversity, absent_post_mapping_summ_stats = calculate_metrics(absent_dataframes_normalized, absent_dataframes)
                st.info("Metrics calculated.")

                # Step 8.2: Combine pre- and postMapping information of metrices, where needed
                combined_metrics = combine_metrics(pre_mapping_metrics, present_post_mapping_metrics, df_type="metrics", dfvalues_are_rel_abundances=dfvalues_are_rel_abundances)
                combined_summ_stats = combine_metrics(pre_mapping_summ_stats, present_post_mapping_summ_stats, df_type="summ_stats", dfvalues_are_rel_abundances=dfvalues_are_rel_abundances)
                st.info("Metrics combined.")

                # Optional Step: If stratification groups are provided, stratify dataframe on these groups & calculate metrices on them too
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
                    st.info(f'Stratifying taxa dataframes using following groups: {stratification_groupnames}.')

                # Step 9: Store all result dataframes in a structure & save, if output-path is provided
                dataframe_groups = {'normalized_preMapped': preMapped_dataframes_normalized,
                                    'renormalized_mapped_forModelling': present_dataframes_adjForModelling,
                                    'normalized_mapped': present_dataframes_normalized,
                                    'normalized_unmapped': absent_dataframes_normalized,
                                    'metrics': [combined_metrics, combined_summ_stats, \
                                                pre_mapping_abundance_metrics, \
                                                present_post_mapping_abundance_metrics, \
                                                absent_post_mapping_abundance_metrics, \
                                                pre_mapping_beta_diversity, \
                                                present_post_mapping_beta_diversity, \
                                                absent_post_mapping_beta_diversity]}

                dataframe_groups.update(stratification_groups)

                # --- Output and Downloads ---
                if output_path is not None:
                    # Create the output directory if it doesn't exist
                    if not os.path.exists(output_path):
                        os.makedirs(output_path)
                    save_dataframes(dataframe_groups, output_path, output_format)

                    renamed_dataframe.to_csv(os.path.join(output_path, 'preprocessedInput_afterRenaming.csv'), sep=',')
                    st.success(f"Output saved to {output_path} in {output_format} format.")

                else:
                    st.warning("No output path specified.  Displaying dataframes for download.")

                    for name, data in dataframe_groups.items():
                        st.subheader(f"{name}")
                        if isinstance(data, list):  #metrics are stored in lists
                            for df in data:
                                st.dataframe(df)
                                st.download_button(
                                    label=f"Download {name} as {output_format}",
                                    data=convert_df(df, output_format),
                                    file_name=f"{name}.{output_format}",
                                    mime=f"text/{output_format}" if output_format != "xlsx" else "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                                )
                        else:
                            st.dataframe(data)
                            st.download_button(
                                label=f"Download {name} as {output_format}",
                                data=convert_df(data, output_format),
                                file_name=f"{name}.{output_format}",
                                mime=f"text/{output_format}" if output_format != "xlsx" else "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            )

            st.success("MARS analysis complete!")

        except Exception as e:
            st.error(f"An error occurred: {e}")
