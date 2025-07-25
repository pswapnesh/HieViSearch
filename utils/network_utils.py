import networkx as nx
import numpy as np
import pandas as pd


def make_network(hdb, df, wt_nan=1e12,min_lambda = -1):
    """
    Constructs a directed graph from HDBSCAN's condensed tree.

    Parameters:
        hdb: HDBSCAN object with condensed_tree_ attribute.
        df (pd.DataFrame): Dataframe containing node attributes.
        wt_nan (float): Weight assigned when lambda_val is NaN.

    Returns:
        nx.DiGraph: A directed graph with nodes and edges from the condensed tree.
    """
    G = nx.DiGraph()

    # Add all nodes from the condensed tree
    all_nodes = set(hdb.condensed_tree_._raw_tree['parent']).union(set(hdb.condensed_tree_._raw_tree['child']))
    G.add_nodes_from(all_nodes)

    # Add edges with weights
    for row in hdb.condensed_tree_._raw_tree:
        parent, child, lambda_val = int(row['parent']), int(row['child']), row['lambda_val']
        weight = lambda_val if np.isfinite(lambda_val) else wt_nan
        if min_lambda > 0:
            if lambda_val > min_lambda:                
                G.add_edge(parent, child, weight=weight)#, distance=1 / weight if weight != 0 else wt_nan)
        else:
            G.add_edge(parent, child, weight=weight)#, distance=1 / weight if weight != 0 else wt_nan)

    # Assign attributes to existing nodes
    for i, (idx, row) in enumerate(df.iterrows()):
        node_id = int(i)
        if node_id in G.nodes:
            G.nodes[node_id].update(row.to_dict())

    return G


import networkx as nx

def get_subtree_from_ancestor(graph: nx.DiGraph, accession: str, number_of_predecessors: int=2) -> nx.DiGraph:
    """
    Given a graph and an accession, find the subtree rooted at the 2nd ancestor of the node with that accession.

    Parameters:
    - graph: A directed NetworkX graph (assumed to be a tree)
    - accession: The accession string to look for in node attributes

    Returns:
    - A subgraph rooted at the 2nd ancestor of the node with the given accession
    """
    # Step 1: Find node ID by accession
    target_node = None
    for node_id, attrs in graph.nodes(data=True):
        if attrs.get("accession") == accession:
            target_node = node_id
            break
    if target_node is None:
        raise ValueError(f"Accession '{accession}' not found in graph.")

    # Step 2: Go to the 2nd ancestor
    current = target_node
    for i in range(2):
        predecessors = list(graph.predecessors(current))        
        if not predecessors:
            #raise ValueError(f"Node '{current}' does not have enough ancestors.")
            break
        current = predecessors[0]  # Assumes single parent (tree)

    ancestor_node = current

    # Step 3: Get all descendants of this ancestor
    descendants = nx.descendants(graph, ancestor_node)
    descendants.add(ancestor_node)  # Include the ancestor itself

    # Step 4: Return the induced subgraph
    return graph.subgraph(descendants).copy()

import networkx as nx

def graph_to_newick(graph: nx.DiGraph) -> str:
    """
    Converts a directed tree (DiGraph) into Newick format using the 'accession'
    attribute as the node label, and writes it to a file.

    Parameters:
    - graph: A directed NetworkX tree
    - output_file: Path to the output Newick file
    """
    # Step 1: Identify the root node (no incoming edges)
    roots = [n for n in graph.nodes if graph.in_degree(n) == 0]
    if len(roots) != 1:
        raise ValueError(f"Graph must have exactly one root, found {len(roots)}.")
    root = roots[0]

    # Step 2: Recursive function to convert to Newick
    def to_newick(node):
        children = list(graph.successors(node))
        label = graph.nodes[node].get("accession", str(node))
        if children:
            return f"({','.join(to_newick(child) for child in children)}){label}"
        else:
            return label

    # Step 3: Convert and write to file
    newick_str = to_newick(root) + ";"
    return newick_str


def convert_similarity_to_distance(graph: nx.Graph, similarity_attr: str = "weight") -> nx.Graph:
    """
    Returns a copy of the graph where edge weights are inverted: distance = 1 / similarity.
    Assumes similarity > 0 for all edges.
    """
    g = graph.copy()
    for u, v, data in g.edges(data=True):
        sim = data.get(similarity_attr, 1.0)
        data["distance"] = 1.0 / sim if sim > 0 else float("inf")
    return g


def merge_graphs_by_accession_safe(graph_list, accession_key= "accession"):
    """
    Merges graphs by the 'accession' attribute. If a node lacks 'accession',
    it is assigned a fallback ID: f"{node_id}_{graph_index}" to ensure uniqueness.

    Parameters:
    - graph_list: List of NetworkX graphs (Graph or DiGraph)
    - accession_key: Node attribute to use as unique identity (default: 'accession')

    Returns:
    - A merged graph where nodes are keyed by accession or fallback ID.
    """
    is_directed = any(isinstance(g, nx.DiGraph) for g in graph_list)
    merged = nx.DiGraph() if is_directed else nx.Graph()

    for i, g in enumerate(graph_list):
        node_id_to_acc = {}

        for node_id, attrs in g.nodes(data=True):
            acc = attrs.get(accession_key)
            if acc is None:
                acc = f"{node_id}"  # fallback unique ID
                attrs[accession_key] = acc  # insert into attributes for consistency
            node_id_to_acc[node_id] = acc

            if merged.has_node(acc):
                merged.nodes[acc].update(attrs)
            else:
                merged.add_node(acc, **attrs)

        for u, v, edge_attrs in g.edges(data=True):
            acc_u = node_id_to_acc.get(u)
            acc_v = node_id_to_acc.get(v)
            if acc_u is None or acc_v is None:
                continue  # should not happen, but guard anyway
            if merged.has_edge(acc_u, acc_v):
                merged[acc_u][acc_v].update(edge_attrs)
            else:
                merged.add_edge(acc_u, acc_v, **edge_attrs)

    return merged

def merge_graphs(graph_list):
    """
    Merges a list of NetworkX graphs into a single graph.
    Nodes with the same ID are merged; attributes like 'accession' are preserved.

    Parameters:
    - graph_list: List of NetworkX graphs (DiGraph or Graph)

    Returns:
    - A new NetworkX graph containing all nodes and edges from the input graphs.
    """
    # Detect if the input is directed
    is_directed = any(isinstance(g, nx.DiGraph) for g in graph_list)
    combined = nx.DiGraph() if is_directed else nx.Graph()

    for g in graph_list:
        for node, attrs in g.nodes(data=True):
            if combined.has_node(node):
                combined.nodes[node].update(attrs)  # merge attributes
            else:
                combined.add_node(node, **attrs)
        
        for u, v, edge_attrs in g.edges(data=True):
            if combined.has_edge(u, v):
                combined[u][v].update(edge_attrs)  # merge edge attributes
            else:
                combined.add_edge(u, v, **edge_attrs)

    return combined


def get_subgraph_for_accessions(G, accessions):
    """
    Extracts the minimal subtree containing the found accessions.
    Silently skips any accessions not found in the graph.

    Parameters:
        G (nx.DiGraph): Directed tree where some leaves have 'accession' attributes.
        accessions (list): List of accession strings to include.

    Returns:
        nx.DiGraph: Subgraph connecting all found accessions to root.
    """
    # Map accession → node
    acc_to_node = {
        data["accession"]: n
        for n, data in G.nodes(data=True)
        if "accession" in data
    }

    # Keep only accessions that are found
    nodes = [acc_to_node[acc] for acc in accessions if acc in acc_to_node]
    if not nodes:
        return G.subgraph([]).copy()  # Return empty graph

    # Collect all paths from each node to root
    nodes_in_subtree = set()
    for node in nodes:
        while True:
            nodes_in_subtree.add(node)
            preds = list(G.predecessors(node))
            if not preds:
                break  # reached root
            node = preds[0]

    return G.subgraph(nodes_in_subtree).copy()

def get_successors(G,isleaf,predecessor):
    if isleaf.get(predecessor):
        return [predecessor]
    successors = list(G.successors(predecessor))    
    #leaves = [isleaf.get(n) for n in successors]
    return successors

def get_local_nodes(G,isleaf,source_node,depth=2):

    current = source_node
    for i in range(depth):        
        predecessors = list(G.predecessors(current))
        if not predecessors:
            #raise ValueError(f"Node '{current}' does not have enough ancestors.")
            break
        predecessor = predecessors[0]
        current = predecessor
    
    successors = get_successors(G,isleaf,current)
    all_nodes = [current] + list(np.ravel(np.array(successors)))
    for i in range(2*depth):
        tmp = [a for s in successors for a in get_successors(G,isleaf,s)] 
        all_nodes += tmp
        successors = tmp
    all_nodes = np.unique(np.array(all_nodes))
    return G.subgraph(all_nodes).copy()