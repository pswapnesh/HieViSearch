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


from utils.prodigal import *
def main(args):
    """
    HieVi main function.
    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    print("Predicting genes... ")
    embeddings = process_multifasta(args.query_fasta_path,args.output_folder,esm_model)

    # open the query vectors    
    query_phage_ids = [emb["accessions"] for emb in embeddings]
    query_mprs = [emb["embedding"] for emb in embeddings]

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
    first_nearest_accessions = [[metadata(i) for i in indices] for indices in first_nearest_indices]    
    
    nearest_df = pd.concat([pd.DataFrame({"query_accession":[d['accession']]*knn,"nearest_accession":first_nearest_accessions[i]}) for i,d in enumerate(embeddings)])
    nearest_df.to_csv(os.path.join(args.output_folder,"HieVi_nearest_accessions.csv"))
    
    first_nearest_means = np.array([np.mean(np.array([index.reconstruct(int(idx)) for idx in indices]),axis = 0) for indices in first_nearest_indices])
    distances, indices = index.search(first_nearest_means, 256)
    unique_indices = np.unique(np.ravel(indices))
    

    emb_db = np.array([index.reconstruct(int(idx)) for idx in unique_indices])
    acc_db = [metadata(i) for i in unique_indices]
    
    all_emb = np.concatenate((emb_db,query_mprs),axis = 0)

    acc_df = pd.DataFrame("accession":acc_db+list(query_phage_ids))
    distances = euclidean_distances(all_emb).astype('float')
    clusterer = hdbscan.HDBSCAN(min_cluster_size = 2,n_jobs = -1,min_samples = 1,allow_single_cluster = False,cluster_selection_method = "leaf",metric = 'precomputed',gen_min_span_tree=True)
    clusterer.fit(distances)

    G = make_network(clusterer,acc_df,min_lambda = -1)
    nx.write_gexf(G,os.path.join(args.output_folder,'network.gexf'))

    # save nearest neighbours in a single pandas. [query_accession, nearest_accession, distance]


    # for the density graph 



    indices = np.unique(np.ravel(indices))
    print(f"Nearest neighbor search completed. Found {len(indices)} unique neighbors.")
    



    

