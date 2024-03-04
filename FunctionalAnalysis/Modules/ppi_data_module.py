import requests

def retrieve_ppi_data(gene_symbols, species="9606", confidence_score=0.7):

    """
    Retrieve Protein-Protein Interaction (PPI) data from the STRING database.

    Parameters:
    - gene_symbols (list): List of gene symbols.
    - species (str): NCBI Taxonomy ID for the species (default: "9606" for Homo sapiens).
    - confidence_score (float): Minimum confidence score for interactions (default: 0.7).

    Returns:
    - ppi_data (list): List of dictionaries containing PPI information.
    """
    
    base_url = "https://string-db.org/api"
    method = "interaction_partners"

    ppi_data = []

    for gene_symbol in gene_symbols:
        try:
            # Make a request to the STRING API
            response = requests.get(
                f"{base_url}/{method}?identifiers={gene_symbol}&species={species}&required_score={confidence_score}"
            )

            if response.status_code == 200:
                # Parse the JSON response
                data = response.json()

                # Extract relevant information from the response
                if data.get("error"):
                    print(f"Error retrieving PPI data for {gene_symbol}: {data['error']}")
                else:
                    ppi_data.append({
                        "Gene Symbol": gene_symbol,
                        "PPI Partners": [partner["preferredName"] for partner in data["partners"]],
                        "Confidence Score": confidence_score,
                    })
            else:
                print(f"Error accessing STRING API for {gene_symbol}: {response.status_code}")
        except Exception as e:
            print(f"Error retrieving PPI data for {gene_symbol}: {str(e)}")

    return ppi_data
