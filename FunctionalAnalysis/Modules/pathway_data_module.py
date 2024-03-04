from bioservices import KEGG

def retrieve_pathway_data(gene_symbols, organism_code="hsa"):

    """
    Retrieve pathway data from the KEGG database.

    Parameters:
    - gene_symbols (list): List of gene symbols.
    - organism_code (str): KEGG organism code (default: "hsa" for Homo sapiens).

    Returns:
    - pathway_data (list): List of dictionaries containing pathway information.
    """
    
    kegg = KEGG()

    pathway_data = []

    for gene_symbol in gene_symbols:
        try:
            # Search for pathways associated with the gene in KEGG
            pathways = kegg.get_pathway_by_gene(gene_symbol, organism=organism_code)

            # Extract relevant information from the pathway data
            for pathway_id, pathway_name in pathways.items():
                pathway_data.append({
                    "Gene Symbol": gene_symbol,
                    "Pathway ID": pathway_id,
                    "Pathway Name": pathway_name,
                })

        except Exception as e:
            print(f"Error retrieving pathway data for {gene_symbol}: {str(e)}")

    return pathway_data
