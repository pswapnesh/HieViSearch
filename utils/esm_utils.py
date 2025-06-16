import torch
import numpy as np
import esm
from Bio import SeqIO
from io import StringIO

class EsmEmbedding:
    def __init__(self, model="3b"):
        """
        Initialize the ESM embedding model based on the specified version.

        Parameters:
        - model (str): The model version to load. Options are '650m', '15b', '3b', or others (default: '15b').
        """
        # Load model and define layers based on the version
        model_map = {
            "650m": (esm.pretrained.esm2_t33_650M_UR50D, [33]),
            "3b": (esm.pretrained.esm2_t36_3B_UR50D, [36]),
            "default": (esm.pretrained.esm2_t36_3B_UR50D, [36]),
        }

        self.device = "cuda"
        model_loader, self.layers = model_map.get(model, model_map["default"])
        self.model, self.alphabet = model_loader()

        # Prepare the batch converter and model
        self.batch_converter = self.alphabet.get_batch_converter()
        self.model.eval()  # Disable dropout for deterministic results
        
        print("########## Model loaded.")
        
        torch.cuda.empty_cache()
        if torch.cuda.is_available():
            self.device = "cuda"
            self.model = self.model.cuda()
            print("Transferred model to GPU")
        else:
            self.device = "cpu"
            print("GPU not available. Running on CPU.")

    def predict(self, data):
        """
        Generate embeddings for the provided data.

        Parameters:
        - data: A batch of sequences to process.

        Returns:
        - np.ndarray: A numpy array containing the embedding vector.
        """
        

        # Convert the input sequence to tokens
        batch_labels, batch_strs, batch_tokens = self.batch_converter(data)
        batch_lens = (batch_tokens != self.alphabet.padding_idx).sum(1)

        batch_tokens = batch_tokens.to(device=self.device, non_blocking=True)

        # Run the model inference
        with torch.no_grad():
            results = self.model(
                batch_tokens, repr_layers=self.layers, return_contacts=False
            )

        # Extract embeddings for the last layer (no dictionary return, just the vector)
        layer = self.layers[-1]  # Using only the last layer
        token_rep = results["representations"][layer].to(torch.float64)

        return token_rep[0, 1 : -1].mean(axis=0),token_rep[0, 0]  # This will return a tensor if kept on GPU

    def predict_proteome(self, protein_fasta):
        """
        Predicts embeddings for each protein sequence in the given FASTA string.
        Loads the model once and processes all sequences within a single GPU context.
        
        Args:
            protein_fasta (str): A FASTA string containing protein sequences.
        Returns:
            np.ndarray: The mean embedding vector across all protein sequences.
        """       
        protein_sequences = {}
        for record in SeqIO.parse(StringIO(protein_fasta), "fasta"):
            protein_sequences[record.id] = str(record.seq)
    
        embeddings = []
        
        # Process all sequences within this single GPU context
        for prot_id, sequence in protein_sequences.items():
            # Use the helper method directly - no GPU context switch
            mean_embedding,_ = self.predict([('name',sequence)])
            embeddings.append(mean_embedding)
    
        # Stack embeddings and normalize
        embeddings_stack = torch.stack(embeddings, dim=0).to(torch.float64)
        embeddings_stack = torch.nn.functional.normalize(embeddings_stack, p=2, dim=1)
        embeddings_mean = torch.mean(embeddings_stack, dim=0)
    
        return embeddings_mean.cpu().numpy()  # Move to CPU and convert to numpy
