# Hievi

# Download dataset and faiss index

*Nearest neighbour index (Faiss index) *
https://sdrive.cnrs.fr/s/kTJyFtzHNAg4t4f

*Metadata*
https://sdrive.cnrs.fr/s/kBwgyYe6rwjmNgZ


# Installation

**Install prodigal**
* This is required, in the command line prodigal -h  must be working before procceeding 
If linux-64 or osx-64, use conda.
conda install bioconda::prodigal
Else check
https://github.com/hyattpd/Prodigal/wiki/installation


conda env create -f environment.yml


# Usage

hievi --experiment_name test --query_fasta_path test/test_contigs.fasta --faiss_index_path /media/microscopie-lcb/swapnesh/protein/embeddings/phages/14Apr2025/Inphared_14Apr_float64_t36_meanCls_3b.zarr_faiss_index.bin --output_folder /media/microscopie-lcb/swapnesh/protein/embeddings/phages/CutTest


