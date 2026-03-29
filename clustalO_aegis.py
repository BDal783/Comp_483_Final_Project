#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 28 2026 by leahbriscoe830

Test of the AEGIS pipeline using the ClustalO alignment tool
"""

import subprocess
import argparse
import sys
import time

def main():

    # Set up parser
    parser = argparse.ArgumentParser(description='Runs the AEGIS pipeline')
    parser.add_argument('-i', '--input', help='input file', required=True)

    # Define infile
    args = parser.parse_args(sys.argv[1:])
    infile = args.input

    # Retrieve data from ncbi_extract.py
    subprocess.run(['python', 'ncbi_extract.py', '-i', infile, '-o', 'example_fasta.txt'], check=True)

    # Time ClustalO alignment with 1000 sequences
    start_time = time.perf_counter()

    # Fix fasta headers for ClustalO input
    with open('example_fasta.txt', 'r') as infile, open('proteinIDs.txt', 'w') as outfile:
    # Iterate through the fasta
        for line in infile:
            # Only modify headers
            if line.startswith(">"):
                # Extract protein_id
                protein_id = line.split(' ')[0]
                outfile.write(str(protein_id)+"\n")
            # Write the sequence to the modified fasta as normal
            else:
                outfile.write(line) 

    with open('clustal_aligned.txt', 'w') as o:
        # Perform MSA using ClustalO, output in fasta format
        subprocess.run(['clustalo', '-i', 'proteinIDs.txt', '-o', o.name, '--outfmt=fasta', '--verbose', '--force', '--quiet'], check=True)

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print("ClustalO Time: " + str(elapsed_time))

if __name__ == '__main__':
    main()