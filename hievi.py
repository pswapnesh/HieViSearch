import zarr
import networkx as nx
import argparse
import os
import faiss 
import pandas as pd
from utils.proteome_process import *
from utils.fasta_utils import *
from utils.esm_utils import *
from utils.network_utils import *
import json
from sklearn.preprocessing import normalize
from sklearn.metrics.pairwise import euclidean_distances, cosine_distances
import hdbscan 
from utils.prodigal import *

esm_model = EsmEmbedding()

def get_argument_parser():
    """
    Create and return the argument parser.
    """
    parser = argparse.ArgumentParser(description="Process and cluster proteome representations.")
    parser.add_argument("--experiment_name", type=str, required=True, help="Experiment name prefix.")
    parser.add_argument("--query_fasta_path", type=str, required=True, help="Experiment name prefix.")
    parser.add_argument("--faiss_index_path", type=str, required=True, help="Path to FAISS index.")
    parser.add_argument("--output_folder", type=str, required=True, help="Output folder path.")
    parser.add_argument("--k_neighbours", type=int, default=16, help="Number of neighbors for FAISS search.")    
    return parser


def main(args):
    """
    HieVi main function.
    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    print("Predicting genes... ")
    embeddings = process_multifasta(args.query_fasta_path,args.output_folder,esm_model)

    # open the query vectors    
    query_phage_ids = [emb["accession"] for emb in embeddings]
    query_mprs = np.array([emb["embedding"] for emb in embeddings])

    number_of_queries = len(query_mprs)

    # if number > 32, and minimum distance between queries is less than some threshold then group them and search may be?

    # 1. Search one by one
    # 2. group and search for each group
    
    # Open and load the JSON file
    json_file_path = args.faiss_index_path.replace("faiss_index.bin","faiss_metadata.json")
    with open(json_file_path, 'r') as f:
        metadata = json.load(f)

    index = faiss.read_index(args.faiss_index_path)
    knn = args.k_neighbours
    first_nearest_distances, first_nearest_indices = index.search(query_mprs, knn)        
    first_nearest_accessions = [[metadata[str(i)] for i in indices] for indices in first_nearest_indices]    
    
    nearest_df = pd.concat([pd.DataFrame({"query_accession":[d['accession']]*knn,"nearest_accession":first_nearest_accessions[i]}) for i,d in enumerate(embeddings)])
    nearest_df.to_csv(os.path.join(args.output_folder,"nearest_neighbour_accessions.csv"))
    
    first_nearest_means = np.array([np.mean(np.array([index.reconstruct(int(idx)) for idx in indices]),axis = 0) for indices in first_nearest_indices])
    distances, indices = index.search(first_nearest_means, 256)
    indices = list(np.ravel(first_nearest_indices)) + list(np.ravel(indices))
    unique_indices = np.unique(np.array(indices))
    print(f"Nearest neighbor search completed. Found {len(unique_indices)} unique neighbors.")
    

    emb_db = np.array([index.reconstruct(int(idx)) for idx in unique_indices])
    acc_db = [metadata[str(i)]  for i in unique_indices]
    
    all_emb = np.concatenate((emb_db,query_mprs),axis = 0)

    
    print("Running density clustering")
    distances = euclidean_distances(all_emb).astype('float')
    clusterer = hdbscan.HDBSCAN(min_cluster_size = 2,n_jobs = -1,min_samples = 1,allow_single_cluster = False,cluster_selection_method = "leaf",metric = 'precomputed',gen_min_span_tree=True)
    clusterer.fit(distances)

    acc_df = pd.DataFrame({"accession":acc_db+list(query_phage_ids),"clust_id": clusterer.labels_})

    print("Making tree")
    G = make_network(clusterer,acc_df,min_lambda = -1)
    G = convert_similarity_to_distance(G)

    # Accession → node
    accession_to_node = {
        data["accession"]: node for node, data in G.nodes(data=True) if "accession" in data
    }
    isleaf = {node: True if G.out_degree(node) == 0 else False for node in G.nodes()}


    nearest_in_tree_df = []
    graphs_list = []    
    nx.write_gexf(G,os.path.join(args.output_folder,'full.gexf'))
    for q,accs in zip(query_phage_ids,first_nearest_accessions):
        subtree = get_local_nodes(G,isleaf,accession_to_node[q],depth = 2)
        nearest_in_tree = pd.DataFrame({'query_accession':[q]*len(accs),"accession":accs})
        graphs_list.append(subtree)
        #nx.write_gexf(subtree,os.path.join(args.output_folder,q,'network.gexf'))        
        nearest_in_tree_df.append(nearest_in_tree)

    nearest_in_tree_df = pd.concat(nearest_in_tree_df)        
    G_full = merge_graphs(graphs_list) #merge_graphs_by_accession_safe(graphs_list)
    nx.write_gexf(G_full,os.path.join(args.output_folder,'hievi_network.gexf'))
    nearest_in_tree_df.to_csv(os.path.join(args.output_folder,args.experiment_name+'_nearest_in_tree.csv'))
    pd.DataFrame({"accession":acc_db}).to_csv(os.path.join(args.output_folder,'relevant_accessions.csv'))


if __name__ == "__main__":
    parser = get_argument_parser()
    args = parser.parse_args()  # NOT sys.argv[1:] manually
    main(args)


