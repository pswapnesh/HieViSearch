# Hievi

Single 

https://huggingface.co/spaces/pswap/hievi

---


## 📥 Download dataset and faiss index

Download and save the two files in the same directory. 

*Nearest neighbour index (Faiss index)*
https://sdrive.cnrs.fr/s/kTJyFtzHNAg4t4f

*Metadata*
https://sdrive.cnrs.fr/s/kBwgyYe6rwjmNgZ



## 🛠️ Installation

### 1. Create and activate the Conda environment

```bash
conda env create -f environment.yml
conda activate hievi
```

Clone this project onto a local folder
```bash
cd your_folder
git clone https://github.com/pswapnesh/HieViSearch.git
```


## 🚀 Usage
```bash
python hievi.py \
  --experiment_name test \
  --query_fasta_path test/test_contigs.fasta \
  --faiss_index_path path/to/faissbin.bin \
  --output_folder path/to/output
```
## Expected execution time
Apart from the model loading, which takes a few seconds. Each contig is processed in less 2 seconds per contig (of size ~ 30kbp) (GPU).

## 📂 Output
The output is structured as bgiven below.
- ```nearest_neighbour_accessions.csv``` contains a table of the query accession and the k-nearest neighbour accessions in the database
- ```experimentName_nearest_in_tree.csv``` contains the k-nearest accession in the condensed tree
- ```hievi_network.gexf``` contains the combined tree of all the qery and relevant accessions in the database 

```php-template
Output_Folder/
├── Phage_01/
├── Phage_02/
│   └── prodigal_outputs
├── hievi_network.gexf
├── <experiment_name>_nearest_in_tree.csv
└── nearest_neighbour_accessions.csv
```

## Visualization
Cytoscape is recommended with gexf-app and yFiles Layout Algorithms.

## 📫 Cite
Swapnesh Panigrahi, Mireille Ansaldi, Nicolas Ginet
https://www.biorxiv.org/content/10.1101/2024.12.17.627486v1