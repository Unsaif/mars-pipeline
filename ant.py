from requests_html import HTMLSession
from tax_id import get_tax_id
import pickle
import streamlit as st


def remove_consecutive_duplicates(lst):
    new_lst = [lst[i] for i in range(len(lst) - 1) if lst[i] != lst[i + 1]]
    new_lst.append(lst[-1])  # Always include the last element
    return new_lst


def homosynonym_finder(species_list):
    # Open a pickle file and load it into a variable
    with open("agora2_species.pkl", "rb") as f:
        agora2_species = pickle.load(f)

    # Initialize an empty dictionary to store species and their Tax IDs
    taxid_dic = {}

    st.info("Finding Taxonomic IDs for given species..")

    progress_text = "Operation in progress. Please wait."
    bar_1 = st.progress(0)

    # For each species in the list, try to get its Tax ID and store it in the dictionary
    for i, species in enumerate(species_list):
        length = len(species_list)
        progress = round(((i + 1) / length) * 100)
        bar_1.progress(progress, text=progress_text)
        try:
            taxid = get_tax_id(species)
            taxid_dic[species] = taxid
            progress_text = f"Found Taxonomic ID for {species}"
        except Exception as e:
            # If unable to find Tax ID for a species, print an error message
            progress_text = f"Unable to find Taxonomic ID for {species}"

    st.success("Done!")

    # Start a session
    session = HTMLSession()

    st.info("Finding potential homosynonyms and mapping to AGORA2 resource..")

    progress_text = "Operation in progress. Please wait."
    bar_2 = st.progress(0)

    # Initialize lists to store species that are already in Agora2 and homosynonyms found in Agora2
    already_in_agora2 = []
    homosynonym_in_agora2 = {}

    # For each species and its corresponding Tax ID in the dictionary
    for i, (species, taxid) in enumerate(taxid_dic.items()):
        length = len(taxid_dic)
        progress = round(((i + 1) / length) * 100)
        bar_2.progress(progress, text=f"Processing {species}..")
        try:
            # Send a GET request to the NCBI Taxonomy Browser
            r = session.get(
                f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={taxid}"
            )

            # Use CSS selector to find the homosynonyms
            homosynonym_elements = r.html.find("div strong i")
            homosynonym_elements = [element.text for element in homosynonym_elements]

            try:
                # Try to remove consecutive duplicates from the list
                homosynonym_elements = remove_consecutive_duplicates(
                    homosynonym_elements
                )
            except:
                # If unable to do so, ignore and move on
                pass

            # If there are any elements left in the list
            if homosynonym_elements:
                # Pair up elements and join each pair together with a space
                potential_homosynonyms = [
                    " ".join(
                        [
                            homosynonym_elements[i],
                            homosynonym_elements[i + 1],
                        ]
                    )
                    for i in range(0, len(homosynonym_elements) - 1, 2)
                ]

                # Initialize a new list to store refined homosynonyms
                refined_potential_homosynonyms = []
                for homosynonym in potential_homosynonyms:
                    if homosynonym != species.strip():
                        # Only add to the list if the homosynonym is not equal to the original species
                        refined_potential_homosynonyms.append(homosynonym)

            # Remove duplicates from the list of refined potential homosynonyms
            refined_potential_homosynonyms = list(set(refined_potential_homosynonyms))

            # Check if the species is not in Agora2
            if species not in agora2_species:
                if refined_potential_homosynonyms:
                    # For each homosynonym, check if it's in Agora2
                    for homosynonym in refined_potential_homosynonyms:
                        if homosynonym in agora2_species:
                            # If it is, add to the dictionary of homosynonyms in Agora2
                            homosynonym_in_agora2[species] = homosynonym
            else:
                # If the species is already in Agora2, add it to the corresponding list
                already_in_agora2.append(species)
        except Exception as e:
            print(e)
            continue

        # Return the list and dictionary
    return already_in_agora2, homosynonym_in_agora2
