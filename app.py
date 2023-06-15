# Import the necessary libraries
import streamlit as st
import pandas as pd
import base64
from ant import homosynonym_finder


def convert_df(df, file_format):
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
        elif "excel" in uploaded_file.type:
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
st.title("ANT - Automated NCBI Taxonomy browser")

# Let the user upload a file
uploaded_file = st.file_uploader(
    "Choose a TXT, CSV or Excel file", type=["txt", "csv", "xlsx"]
)

# Convert the uploaded file to a list
species_list = file_to_list(uploaded_file)

# If the species list is not empty, find homosynonyms
if species_list is not None:
    species_in_agora2, homosynonyms = homosynonym_finder(species_list)
    st.success("Done!")

    species_df = pd.DataFrame(species_in_agora2, columns=["Species in AGORA2"])
    homosynonyms_df = pd.DataFrame(
        list(homosynonyms.items()),
        columns=["Original Name", "Homosynonym in AGORA2"],
    )
    st.dataframe(species_df)
    st.dataframe(homosynonyms_df)

    # Let the user choose the file format
    file_format = st.selectbox("Select file format", ["csv", "xlsx"])

    # Convert the DataFrame to base64
    b64_species = convert_df(species_df, file_format)
    b64_homosynonyms = convert_df(homosynonyms_df, file_format)

    # Create two columns
    col1, col2 = st.columns(2)

    # Create the download button
    with col1:
        st.download_button(
            label="Download Species",
            data=b64_species,
            file_name=f"data.{file_format}",
            mime=f"text/{file_format}",
        )

    with col2:
        # Create the download button
        st.download_button(
            label="Download Homosynonyms",
            data=b64_homosynonyms,
            file_name=f"data.{file_format}",
            mime=f"text/{file_format}",
        )
