# Hievi: 

Single 

https://huggingface.co/spaces/pswap/hievi

---


##  Download dataset and faiss index

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

## 🚀 Usage
```bash
python hievi.py \
  --experiment_name test \
  --query_fasta_path test/test_contigs.fasta \
  --faiss_index_path path/to/faissbin.bin \
  --output_folder path/to/output
```

## 📂 Output

## 📫 Cite