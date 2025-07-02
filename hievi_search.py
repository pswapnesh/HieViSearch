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
from scipy.spatial.distance import cdist
from utils.plotter import *

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
    
    # Open and load the JSON file
    json_file_path = args.faiss_index_path.replace("faiss_index.bin","faiss_metadata.json")
    with open(json_file_path, 'r') as f:
        metadata = json.load(f)

    index = faiss.read_index(args.faiss_index_path)
    knn = args.k_neighbours
    
    print("Searching ... ")
    first_nearest_distances, first_nearest_indices = index.search(query_mprs, knn)        
    first_nearest_accessions = [[metadata[str(i)] for i in indices] for indices in first_nearest_indices]    
    
    nearest_df = pd.concat([pd.DataFrame({"query_accession":[d['accession']]*knn,"nearest_accession":first_nearest_accessions[i],"distance":first_nearest_distances[i]}) for i,d in enumerate(embeddings)])
    nearest_df.to_csv(os.path.join(args.output_folder,"nearest_neighbour_accessions.csv"))


if __name__ == "__main__":
    parser = get_argument_parser()
    args = parser.parse_args()  # NOT sys.argv[1:] manually
    main(args)


