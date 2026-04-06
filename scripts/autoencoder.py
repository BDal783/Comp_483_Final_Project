
from __future__ import annotations
import numpy as np
import pandas as pd
# from sklearn.model_selection import train_test_split
# from tensorflow.keras.models import Model
# from tensorflow.keras.layers import Input, Dense, Flatten, Reshape, Dropout
# import tensorflow as tf
import matplotlib.pyplot as plt
# from scipy.stats import mode
from Bio import SeqIO
import argparse
import sys


def fasta_to_dataframe(fasta_path: str) -> pd.DataFrame:
    """Turns multi-fasta file into pandas dataframe with rows as samples
    and columns as postions in the sequence. Fasta file must be aligned 
    beforehand.

    :param fasta_path: Fasta file to turn into dataframe
    :return: dataframe
    """

    # Create a list of SeqIO objects for each fasta record
    records = list(SeqIO.parse(fasta_path, "fasta"))

    # For each fasta record take the sequence and place in a list
    # Ex: [['A', 'G', 'P'...], ['A', 'C', 'P'....],....]
    seq_data = [list(str(record.seq)) for record in records]

    # Turn into a dataframe and add the id
    df = pd.DataFrame(seq_data)
    df.insert(0, "ID", [record.id for record in records])

    # Set up column names for the dataframe
    # Ex: ['ID', 'Letter1', 'Letter2', 'Letter3'....] Until max letter
    column_names = ['ID']
    lengths = []
    max_letter = len(seq_data[0])
    for i in range(max_letter):
        column_names.append(f"Letter {i+1}")
    df.columns = column_names
    return df


def main():

    # Set up parser
    parser = argparse.ArgumentParser(description='Runs the autoencoder for the AEGIS pipeline')
    parser.add_argument('-i', '--input', help='input file', required=True)

    # Grab input file
    args = parser.parse_args(sys.argv[1:])
    infile = args.input

    # Turn input file into dataframe
    df = fasta_to_dataframe(infile)

    # Get the amino acid sequences as input data
    sequences = df.iloc[:, 1:].astype(str).values

    # Find unique symbols in amino acids
    # Should be 22 symbols with X and '-' included
    amino_acids = np.unique(sequences)

    # Map symbols to integers
    amino_acid_to_int = {aa: i for i, aa in enumerate(amino_acids)}

    # Find length of sequences (Assumes all to be the same due to MSA) and number of symbols 
    input_dim = sequences.shape[1]
    num_amino_acids = len(amino_acids)

    # Initialize a 3D array of zeros: (num_sequences, sequence_length, num_amino_acids)
    # This will hold the one-hot encoded sequences
    encoded_sequences = np.zeros((len(sequences), input_dim, num_amino_acids), dtype=int)

    # Loop over sequences and positions to fill in one-hot vectors
    for i, seq in enumerate(sequences):       # i = sequence index (row)
        for j, aa in enumerate(seq):          # j = position index in the sequence (col)
            # Set the corresponding amino acid index to 1 (encoding vector)
            encoded_sequences[i, j, amino_acid_to_int[aa]] = 1

if __name__ == '__main__':
    main()
