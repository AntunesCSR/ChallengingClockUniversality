from Bio import Entrez
import json


def retrieve_gene_annotation(gene_symbols, email):
    """
    Retrieve gene annotation using BioPython and NCBI Entrez.

    Parameters:
    - gene_symbols (list): List of gene symbols.
    - email (str): Your email address (required for NCBI Entrez).

    Returns:
    - gene_annotation_data (list): List of dictionaries containing gene annotation information.
    """
    Entrez.email = 'antunescsr@gmail.com'

    gene_annotation_data = []

    for gene_symbol in gene_symbols:
        search_results = {}  # Initialize before the try block
        try:
            # Search for gene information
            search_handle = Entrez.esearch(db="gene", term=f"{gene_symbol}[Gene Symbol]")
            search_results = Entrez.read(search_handle)

            esearch_result = search_results.get('esearchresult', {})
            if esearch_result.get('count', 0) > 0:
                gene_id = esearch_result['idlist'][0]

                # Fetch gene information
                gene_handle = Entrez.esummary(db="gene", id=gene_id, retmode="xml")
                gene_info = Entrez.parse(gene_handle)

                # Extract information from the generator
                gene_info_dict = next(gene_info)

                gene_annotation_data.append({
                    "Gene Symbol": gene_symbol,
                    "Gene ID": gene_id,
                    "Description": gene_info_dict['result'][gene_id]['description'],
                    "Biotype": gene_info_dict['result'][gene_id]['type'],
                    "Chromosome": gene_info_dict['result'][gene_id]['chromosome'],
                    # Add more fields as needed
                })
            else:
                print(f"No matching gene found for {gene_symbol}")
        except Exception as e:
            print(f"Error retrieving gene annotation for {gene_symbol}: {str(e)}")
            # Print or log the raw response for further inspection
            print("Raw response:", search_results)

    return gene_annotation_data
