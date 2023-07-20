# Import the necessary libraries
import streamlit as st
import pandas as pd
from ant import ant
from io import BytesIO


def convert_df(df, file_format):
    if file_format == "xlsx":
        output = BytesIO()
        with pd.ExcelWriter(output, engine="openpyxl") as writer:
            df.to_excel(writer, index=False)
        data = output.getvalue()
        return data
    elif file_format == "txt":
        return df.to_csv(index=False, sep="\t").encode("utf-8")
    else:
        return df.to_csv(index=False).encode("utf-8")


def file_to_list(uploaded_file):
    # Reads a file and returns the data as a list
    if uploaded_file is not None:
        file_details = {
            "FileName": uploaded_file.name,
            "FileType": uploaded_file.type,
            "FileSize": uploaded_file.size,
        }
        st.write(file_details)

        if uploaded_file.type == "text/plain":
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


# Add a title to your Streamlit app
st.title(":ant: ANT - Automated NCBI Taxonomy")

st.markdown(
    f"""

##### ANT is a dynamic tool designed to simplify and automate the process of mapping species to a resource selected by the user or the default [AGORA2 resource](https://www.nature.com/articles/s41587-022-01628-0). 

###### By applying reliable web scraping techniques, ANT automates the search for species homosynonyms (alternative taxonomic names) within the [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi).

This app generates comprehensive results, which include:
- All species present in the resource
- Those with recognized homosynonyms in the resource
- The identified NCBI taxonomic IDs
- All discovered homosynonyms

For the user's convenience, these results can be downloaded in your chosen file format, facilitating further exploration and analysis.
"""
)

# Create two columns
col1, col2 = st.columns(2)

# Let the user upload a file
uploaded_file = col1.file_uploader(
    "Choose a TXT, CSV or Excel file ([example file](https://gitlab.com/timhunuig/ant/-/raw/master/example.txt))",
    type=["txt", "csv", "xlsx"],
)
uploaded_resource_file = col2.file_uploader(
    "(Optional) Upload your own resource file",
    type=["txt", "csv", "xlsx"],
)

# Convert the uploaded file to a list
species_list = file_to_list(uploaded_file)
resource = file_to_list(uploaded_resource_file)

st.divider()

# If the species list is not empty, find homosynonyms

if species_list is not None:
    if st.button("Start ANT"):
        species_in_resource, homosynonyms, ncbi_tax_id, all_found_homosynoyms = ant(
            species_list, resource
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

        col1, col2 = st.columns(2)

        col1.dataframe(species_df)
        col2.dataframe(homosynonyms_df, hide_index=True)

        st.divider()

        col1, col2 = st.columns(2)

        col1.dataframe(ncbi_tax_id_df, hide_index=False)
        col2.dataframe(all_found_homosynoyms_df, hide_index=True)

        st.divider()

        # Create two columns
        col1, col2 = st.columns(2)

        # Let the user choose the file format
        file_format = col2.selectbox("Select file format", ["csv", "txt", "xlsx"])

        # Convert the DataFrame to downloadable
        species = convert_df(species_df, file_format)
        homosynonyms = convert_df(homosynonyms_df, file_format)
        ncbi_tax_ids = convert_df(ncbi_tax_id_df, file_format)
        all_found_homosynoyms = convert_df(all_found_homosynoyms_df, file_format)

        # Create the download button
        col1.download_button(
            label="Download Species",
            data=species,
            file_name=f"species.{file_format}",
            mime=f"text/{file_format}",
            disabled=species_df.empty,
        )

        col1.download_button(
            label="Download Homosynonyms",
            data=homosynonyms,
            file_name=f"homosynonyms.{file_format}",
            mime=f"text/{file_format}",
            disabled=homosynonyms_df.empty,
        )

        col1.download_button(
            label="Download NCBI Tax IDs",
            data=ncbi_tax_ids,
            file_name=f"ncbi_tax_ids.{file_format}",
            mime=f"text/{file_format}",
            disabled=ncbi_tax_id_df.empty,
        )

        col1.download_button(
            label="Download All Found Homosynonyms",
            data=all_found_homosynoyms,
            file_name=f"all_found_homosynonyms.{file_format}",
            mime=f"text/{file_format}",
            disabled=all_found_homosynoyms_df.empty,
        )
