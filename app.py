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

st.write(
    f"""

ANT is a tool developed to simplify and automate the process of mapping species to the [AGORA2 resource](https://www.nature.com/articles/s41587-022-01628-0). Leveraging web scraping techniques, 
ANT automates the search for species homosynonyms (alternative taxonomic names) in the [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi).

The results, which include all species present in AGORA2 and those that have homosynonyms in AGORA2, are displayed. 
For convenience, these results can be downloaded in your preferred file format for further use.

"""
)

# Let the user upload a file
uploaded_file = st.file_uploader(
    "Choose a TXT, CSV or Excel file ([example file](https://gitlab.com/timhunuig/ant/-/raw/master/example.txt))",
    type=["txt", "csv", "xlsx"],
)

# Convert the uploaded file to a list
species_list = file_to_list(uploaded_file)

st.divider()

# If the species list is not empty, find homosynonyms
if species_list is not None:
    species_in_agora2, homosynonyms = ant(species_list)
    st.success("Done!")

    st.divider()

    # Create two columns
    col1, col2 = st.columns(2)

    col1.metric("Number of species found in AGORA2", len(species_in_agora2))
    col2.metric("Number of homosynonyms found in AGORA2", len(homosynonyms))

    st.divider()

    species_df = pd.DataFrame(species_in_agora2, columns=["Species in AGORA2"])
    homosynonyms_df = pd.DataFrame(
        list(homosynonyms.items()),
        columns=["Original Name", "Homosynonym in AGORA2"],
    )
    st.dataframe(species_df)
    st.dataframe(homosynonyms_df)

    st.divider()

    # Let the user choose the file format
    file_format = st.selectbox("Select file format", ["csv", "txt", "xlsx"])

    # Convert the DataFrame to base64
    species = convert_df(species_df, file_format)
    homosynonyms = convert_df(homosynonyms_df, file_format)

    # Create two columns
    col1, col2 = st.columns(2)

    # Create the download button
    col1.download_button(
        label="Download Species",
        data=species,
        file_name=f"species.{file_format}",
        mime=f"text/{file_format}",
    )

    col2.download_button(
        label="Download Homosynonyms",
        data=homosynonyms,
        file_name=f"homosynonyms.{file_format}",
        mime=f"text/{file_format}",
    )
