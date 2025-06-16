import networkx as nx
import numpy as np


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
                G.add_edge(parent, child, weight=weight, distance=1 / weight if weight != 0 else wt_nan)
        else:
            G.add_edge(parent, child, weight=weight, distance=1 / weight if weight != 0 else wt_nan)

    # Assign attributes to existing nodes
    for i, (idx, row) in enumerate(df.iterrows()):
        node_id = int(i)
        if node_id in G.nodes:
            G.nodes[node_id].update(row.to_dict())

    return G