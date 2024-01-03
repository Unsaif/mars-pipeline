# Import the necessary libraries
import streamlit as st
import pandas as pd
from ANT.ant import ant
from io import BytesIO

st.set_page_config(layout="wide")

# streamlit_app.py

from MARS.utils import read_file_as_dataframe, merge_files, normalize_dataframes, combine_metrics
from MARS.operations import split_taxonomic_groups, rename_taxa, calculate_metrics, check_presence_in_agora2

def convert_df(df, file_format, index=False):
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
    # Reads a file and returns the data as a list
    if uploaded_file is not None:
        file_details = {
            "FileName": uploaded_file.name,
            "FileType": uploaded_file.type,
            "FileSize": uploaded_file.size,
        }
        st.write(file_details)

        if uploaded_file.type == "text/plain" or uploaded_file.type == "text/tab-separated-values":
            # For txt files, we simply read each line into a list
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
    
st.image("images/logo2.png", use_column_width=True)

st.title("MARS")

st.write("""
""")

col1, col2, col3, col4 = st.columns(4)
uploaded_file1 = col1.file_uploader("Upload Taxonomy Table", type=['csv', 'tsv', 'txt', 'xlsx'])
uploaded_file2 = col2.file_uploader("Upload Feature Table", type=['csv', 'tsv', 'txt', 'xlsx'])
uploaded_file3 = col3.file_uploader("Upload Combined Table", type=['csv', 'tsv', 'txt', 'xlsx'])
uploaded_resource_file = col4.file_uploader(
    "(Optional) Upload your own resource file",
    type=["txt", "csv", "xlsx"],
    )

# Sidebar for optional parameters
st.sidebar.header("Optional Parameters")
cutoff = st.sidebar.slider("Cutoff", 0.0, 100.0, 0.0)
output_format = st.sidebar.radio("Output Format", ["csv", "txt", "xlsx"])
stratification_file = st.sidebar.file_uploader("Stratification File", type=['csv', 'txt', 'xlsx'])
skip_ant = st.sidebar.checkbox('Skip ANT')

if (uploaded_file1 and uploaded_file2) or uploaded_file3:

    # When the button is clicked, the function is executed
    if st.button("Process Files"):
        with st.spinner("Merging files..."):
            if uploaded_file3:
                merged_dataframe = read_file_as_dataframe(uploaded_file3, 0)
                st.success("Files already merged!")
            else:
                merged_dataframe = merge_files(uploaded_file1, uploaded_file2)
                st.success("Merging files: Success!")

        with st.spinner("Splitting taxonomic groups..."):
            taxonomic_dataframes = split_taxonomic_groups(merged_dataframe)
        st.success("Splitting taxonomic groups: Success!")

        with st.spinner("Renaming taxa..."):
            renamed_dataframes = rename_taxa(taxonomic_dataframes)
        st.success("Renaming taxa: Success!")

        # ANT
        if not skip_ant:
            # Convert the uploaded file to a list
            resource = file_to_list(uploaded_resource_file)

            # If the species list is not empty, find homosynonyms

            st.divider()

            species = list(renamed_dataframes['Species'].index)
            species = [sub.replace('_', ' ') for sub in species]

            species_in_resource, homosynonyms, ncbi_tax_id, all_found_homosynoyms = ant(
                species, resource
            )
            st.success("Done!")

            st.divider()

            # Create two columns
            col1, col2 = st.columns(2)

            col1.metric("Number of species found in resource", len(species_in_resource))
            col2.metric("Number of homosynonyms found in resource", len(homosynonyms))

            st.divider()

            species_df = pd.DataFrame(species_in_resource, columns=["Species in Resource"])
            homosynonyms_df = pd.DataFrame(
                list(homosynonyms.items()),
                columns=["Original Name", "Homosynonym in Resource"],
            )
            ncbi_tax_id_df = pd.DataFrame(
                list(ncbi_tax_id.items()), columns=["Species", "NCBI Tax ID"]
            )
            all_found_homosynoyms_df = pd.DataFrame(
                list(all_found_homosynoyms.items()),
                columns=["Species", "All Found Homosynonymns"],
            )

            tab1, tab2, tab3, tab4 = st.tabs(
                [
                    "Species in Resource",
                    "Homosynonyms in Resource",
                    "NCBI Taxonomic IDs",
                    "All Found Homosynonyms",
                ]
            )

            with tab1:
                col1, col2 = st.columns(2)
                col1.dataframe(species_df, use_container_width=True)

                # Convert the DataFrame to downloadable
                species = convert_df(species_df, output_format)

                # Create the download button
                col2.download_button(
                    label="Download Species",
                    data=species,
                    file_name=f"species.{output_format}",
                    mime=f"text/{output_format}",
                    disabled=species_df.empty,
                )
            with tab2:
                col1, col2 = st.columns(2)
                col1.dataframe(
                    homosynonyms_df,
                    hide_index=True,
                    use_container_width=True,
                )

                # Convert the DataFrame to downloadable
                homosynonyms_conv = convert_df(homosynonyms_df, output_format)

                # Create the download button
                col2.download_button(
                    label="Download Homosynonyms",
                    data=homosynonyms_conv,
                    file_name=f"homosynonyms.{output_format}",
                    mime=f"text/{output_format}",
                    disabled=homosynonyms_df.empty,
                )
            with tab3:
                col1, col2 = st.columns(2)
                col1.dataframe(
                    ncbi_tax_id_df,
                    hide_index=False,
                    use_container_width=True,
                )

                # Convert the DataFrame to downloadable
                ncbi_tax_ids = convert_df(ncbi_tax_id_df, output_format)

                # Create the download button
                col2.download_button(
                    label="Download NCBI Tax IDs",
                    data=ncbi_tax_ids,
                    file_name=f"ncbi_tax_ids.{output_format}",
                    mime=f"text/{output_format}",
                    disabled=ncbi_tax_id_df.empty,
                )
            with tab4:
                col1, col2 = st.columns(2)
                col1.dataframe(
                    all_found_homosynoyms_df,
                    hide_index=True,
                    use_container_width=True,
                )

                # Convert the DataFrame to downloadable
                all_found_homosynoyms = convert_df(all_found_homosynoyms_df, output_format)

                # Create the download button
                col2.download_button(
                    label="Download All Found Homosynonyms",
                    data=all_found_homosynoyms,
                    file_name=f"all_found_homosynonyms.{output_format}",
                    mime=f"text/{output_format}",
                    disabled=all_found_homosynoyms_df.empty,
                )

            # rename any new found homosynonyms 
            if homosynonyms:
                with st.spinner("Renaming found homosynonyms..."):
                    homosynonyms = {key.replace(' ', '_'): value.replace(' ', '_') for key, value in homosynonyms.items()}
                    renamed_dataframes["Species"].rename(index=homosynonyms)
                st.success("Renaming found homosynonyms: Success!")

        with st.spinner("Checking presence in AGORA2..."):
            present_dataframes, absent_dataframes = check_presence_in_agora2(renamed_dataframes)
        st.success("Checking presence in AGORA2: Success!")

        with st.spinner("Normalizing dataframes..."):
            normalized_dataframes = normalize_dataframes(renamed_dataframes, cutoff=cutoff)
            normalized_present_dataframes, normalized_absent_dataframes = normalize_dataframes(present_dataframes, cutoff=cutoff), normalize_dataframes(absent_dataframes, cutoff=cutoff)
        st.success("Normalizing dataframes: Success!")

        with st.spinner("Calculating metrics..."):
            pre_agora2_check_metrics = calculate_metrics(renamed_dataframes)
            post_agora2_check_metrics = calculate_metrics(present_dataframes)
            combined_metrics = combine_metrics(pre_agora2_check_metrics, post_agora2_check_metrics)
        st.success("Calculating metrics: Success!")

        stratification_groups = {}
        if stratification_file is not None:
            with st.spinner("Processing stratification groups..."):
                stratification = pd.read_csv(stratification_file)
                for group in stratification["group"].unique():
                    group_columns = list(stratification[stratification["group"] == group]["samples"])
                    pre_group_metrics = calculate_metrics(renamed_dataframes, group=group_columns)
                    post_group_metrics = calculate_metrics(present_dataframes, group=group_columns)
                    combined_group_metrics = combine_metrics(pre_group_metrics, post_group_metrics)
                    group_name = f"{group.lower()}_metrics"
                    stratification_groups[group_name] = combined_group_metrics
            st.success("Processing stratification groups: Success!")

        dataframe_groups = {'normalized': normalized_dataframes, 
                            'present': normalized_present_dataframes, 
                            'absent': normalized_absent_dataframes,
                            'metrics': combined_metrics
                            }

        dataframe_groups.update(stratification_groups)

        endtab1, endtab2, endtab3, endtab4 = st.tabs(
            [
                "Normalized",
                "Present",
                "Absent",
                "Metrics",
            ]
        )
        levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        with endtab1:

            tabs = st.tabs(levels)
            for i, (group, df) in enumerate(dataframe_groups['normalized'].items()):
                
                # Convert the DataFrame to downloadable
                df_conv = convert_df(df, output_format, index=True)

                with tabs[i]:
                    col1, col2 = st.columns(2)
                    col1.dataframe(
                        df,
                        hide_index=False,
                        use_container_width=True,
                    )

                    col2.download_button(
                        label=f"Download {group}",
                        data=df_conv,
                        file_name=f"{group}_normalized.{output_format}",
                        mime=f"text/{output_format}",
                        disabled=df.empty,
                        key=f'{group}_normalized'
                    )
        with endtab2:

            tabs = st.tabs(levels)
            for i, (group, df) in enumerate(dataframe_groups['present'].items()):
                
                # Convert the DataFrame to downloadable
                df_conv = convert_df(df, output_format)

                with tabs[i]:
                    col1, col2 = st.columns(2)
                    col1.dataframe(
                        df,
                        hide_index=False,
                        use_container_width=True,
                    )

                    col2.download_button(
                        label=f"Download {group}",
                        data=df_conv,
                        file_name=f"{group}_present.{output_format}",
                        mime=f"text/{output_format}",
                        disabled=df.empty,
                        key=f'{group}_present'
                    )
        with endtab3:

            tabs = st.tabs(levels)
            for i, (group, df) in enumerate(dataframe_groups['absent'].items()):
                
                # Convert the DataFrame to downloadable
                df_conv = convert_df(df, output_format)

                with tabs[i]:
                    col1, col2 = st.columns(2)
                    col1.dataframe(
                        df,
                        hide_index=False,
                        use_container_width=True,
                    )

                    col2.download_button(
                        label=f"Download {group}",
                        data=df_conv,
                        file_name=f"{group}_absent.{output_format}",
                        mime=f"text/{output_format}",
                        disabled=df.empty,
                        key=f'{group}_absent'
                    )
        with endtab4:

            metric_tabs = st.tabs(['Read Counts', 'Shannon Index', 'Firmicutes/Bacteroidetes Ratio'])

            for i, (group, metrics) in enumerate(dataframe_groups['metrics'].items()):
                for j, (metric, df) in enumerate(metrics.items()):
                    with metric_tabs[j]:
                        metric_tabs[j].subheader(group)
                        # Convert the DataFrame to downloadable
                        df_conv = convert_df(df, output_format)

                        col1, col2 = st.columns(2)
                        col1.dataframe(
                            df,
                            hide_index=False,
                            use_container_width=True,
                        )

                        col2.download_button(
                            label=f"Download {metric}",
                            data=df_conv,
                            file_name=f"{group}_{metric}.{output_format}",
                            mime=f"text/{output_format}",
                            disabled=df.empty,
                            key=f'{group}_{metric}'
                        )