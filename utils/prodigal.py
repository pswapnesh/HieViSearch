import subprocess
import tempfile
import os
from io import StringIO
from Bio import SeqIO
from tqdm import tqdm
import numpy as np

def process_multifasta(multifasta_path,output_folder,esm_model):
    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)
    embeddings = []
    # Parse the multifasta file
    for record in tqdm(SeqIO.parse(multifasta_path, "fasta")):
        seq_name = record.id
        dna_sequence = str(record.seq)

        # Create subfolder for the sequence
        seq_folder = os.path.join(output_folder, seq_name)
        os.makedirs(seq_folder, exist_ok=True)

        # Write individual FASTA file
        fasta_path = os.path.join(seq_folder, f"{seq_name}.fasta")
        fasta_string = f">{seq_name}\n{dna_sequence}\n"
        
        # Predict genes using Prodigal
        cleaned_protein_fasta = predict_genes_prodigal(fasta_string, seq_folder)

        # Predict embeddings
        embeddings_all = esm_model.predict_proteome(cleaned_protein_fasta)
        embedding = np.mean(embeddings_all, axis=0)
        np.save(os.path.join(seq_folder,"embedding.npy"),embedding)
        embeddings.append({"accession": seq_name,"embedding":embedding,"embeddings_all":embeddings_all})
    return embeddings


def predict_genes_prodigal(fasta_string,seq_folder):
    """
    Predicts protein-coding sequences from a given nucleotide FASTA string using Prodigal.
    Args:
        fasta_string (str): A string containing the nucleotide FASTA sequence.
    Returns:
        tuple: A tuple containing:
            - str: The FASTA string of the predicted protein sequences, with stop codons removed.
            - str: The path to the temporary directory.
    """
    # Parse FASTA string to ensure it's valid
    fasta_io = StringIO(fasta_string)
    try:
        record = next(SeqIO.parse(fasta_io, "fasta"))
        sequence = str(record.seq)
        sequence_id = record.id or "sequence"
    except Exception as e:
        return f"Error parsing FASTA string: {str(e)}"

    # Create a temporary directory
    temp_dir = seq_folder
    temp_input_path = os.path.join(temp_dir, "input.fna")
    temp_output_path = os.path.join(temp_dir, "proteins.faa")
    temp_genes_path = os.path.join(temp_dir, "genes.fna")
    temp_annotations_path = os.path.join(temp_dir,"annotations.gbk")

    # Write the input FASTA to a file in the temporary directory
    with open(temp_input_path, 'w') as temp_input:
        temp_input.write(f">{sequence_id}\n{sequence}\n")

    # Run Prodigal via subprocess
    try:
        result = subprocess.run(
            ['prodigal', '-i', temp_input_path, '-a', temp_output_path, '-d', temp_genes_path, '-o', temp_annotations_path, '-p', 'meta'],
            capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        # DO NOT unlink here.  The user wants the temp dir.
        # os.unlink(temp_input_path)
        # os.unlink(temp_output_path)
        return f"Error running Prodigal: {e.stderr}"

    # Read Prodigal output
    try:
        with open(temp_output_path, 'r') as f:
            protein_fasta = f.read()
    except Exception as e:
        # DO NOT unlink here.  The user wants the temp dir.
        # os.unlink(temp_input_path)
        # os.unlink(temp_output_path)
        return f"Error reading Prodigal output: {str(e)}"

    if not protein_fasta.strip():
        return "No proteins predicted." # Return the path even if no proteins

    # Process the output to remove '*' (stop codons)
    cleaned_fasta = []
    fasta_lines = protein_fasta.splitlines()
    for line in fasta_lines:
        if line.startswith(">"):
            cleaned_fasta.append(line)  # Keep header lines
        else:
            # Remove '*' from sequence lines
            cleaned_line = line.replace("*", "")
            if cleaned_line:  # Only append non-empty sequences
                cleaned_fasta.append(cleaned_line)

    cleaned_protein_fasta = "\n".join(cleaned_fasta)

    return cleaned_protein_fasta


