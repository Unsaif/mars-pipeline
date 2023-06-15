from Bio import Entrez


def get_tax_id(species):
    """Get the taxid of a species from ncbi taxonomy by passing the species name to esearch."""
    species = species.replace(" ", "+").strip()
    search = Entrez.esearch(term=species, db="taxonomy", retmode="xml")
    record = Entrez.read(search)
    return record["IdList"][0]


def get_tax_data(taxid):
    """Fetch the taxonomy record of a species using its taxid."""
    search = Entrez.efetch(id=taxid, db="taxonomy", retmode="xml")
    return Entrez.read(search)


# Provide your email to Entrez
Entrez.email = "your_email@example.com"
