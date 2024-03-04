from Bio import Entrez

def retrieve_gene_annotation(gene_symbols):

    """
    Retrieve gene annotation using BioPython and NCBI Entrez.

    Parameters:
    - gene_symbols (list): List of gene symbols.

    Returns:
    - gene_annotation_data (list): List of dictionaries containing gene annotation information.
    """
    
    Entrez.email = "your.email@example.com"  # Provide your email for NCBI Entrez

    gene_annotation_data = []

    for gene_symbol in gene_symbols:
        try:
            # Search for gene information
            search_handle = Entrez.esearch(db="gene", term=f"{gene_symbol}[Gene Symbol]", retmode="json")
            search_results = Entrez.read(search_handle)

            if search_results['esearchresult']['count'] > 0:
                gene_id = search_results['esearchresult']['idlist'][0]

                # Fetch gene information
                gene_handle = Entrez.esummary(db="gene", id=gene_id, retmode="json")
                gene_info = Entrez.read(gene_handle)

                gene_annotation_data.append({
                    "Gene Symbol": gene_symbol,
                    "Gene ID": gene_id,
                    "Description": gene_info['result'][gene_id]['description'],
                    "Biotype": gene_info['result'][gene_id]['type'],
                    "Chromosome": gene_info['result'][gene_id]['chromosome'],
                    # Add more fields as needed
                })
            else:
                print(f"No matching gene found for {gene_symbol}")
        except Exception as e:
            print(f"Error retrieving gene annotation for {gene_symbol}: {str(e)}")

    return gene_annotation_data