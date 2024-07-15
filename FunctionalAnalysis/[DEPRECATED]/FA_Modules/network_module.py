import networkx as nx
import matplotlib.pyplot as plt

def construct_biological_network(gene_annotation_data, ppi_data, pathway_data):

    """
    Construct a biological network using gene annotation, PPI, and pathway data.

    Parameters:
    - gene_annotation_data (list): List of dictionaries containing gene annotation information.
    - ppi_data (list): List of dictionaries containing PPI information.
    - pathway_data (list): List of dictionaries containing pathway information.

    Returns:
    - G (networkx.Graph): Constructed biological network.
    """

    G = nx.Graph()

    # Add nodes from gene annotation data
    for gene_info in gene_annotation_data:
        G.add_node(gene_info["Gene Symbol"], label="Gene")

    # Add edges from PPI data
    for interaction in ppi_data:
        gene_symbol = interaction["Gene Symbol"]
        ppi_partners = interaction["PPI Partners"]

        for partner in ppi_partners:
            G.add_edge(gene_symbol, partner, label="PPI")

    # Add edges from pathway data
    for pathway_info in pathway_data:
        gene_symbol = pathway_info["Gene Symbol"]
        pathway_id = pathway_info["Pathway ID"]

        G.add_edge(gene_symbol, pathway_id, label="Pathway")

    return G




def perform_network_analysis(G):

    """
    Perform basic network analysis on the constructed biological network.

    Parameters:
    - G (networkx.Graph): Constructed biological network.

    Returns:
    - network_analysis_results (dict): Dictionary containing network analysis results.
    """
    
    network_analysis_results = {}

    # Basic network metrics
    network_analysis_results["Number of Nodes"] = G.number_of_nodes()
    network_analysis_results["Number of Edges"] = G.number_of_edges()
    network_analysis_results["Node Degree"] = dict(G.degree())

    # Visualization (optional)
    pos = nx.spring_layout(G)  # Change layout as needed
    nx.draw(G, pos, with_labels=True, font_weight="bold", node_color="skyblue", edge_color="gray", alpha=0.7)
    plt.show()

    return network_analysis_results