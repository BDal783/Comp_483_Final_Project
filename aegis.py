

import subprocess
import argparse
import sys
import time

def main():

    # Set up parser
    parser = argparse.ArgumentParser(description='Runs the AEGIS pipeline')
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-p', '--threshold', help='percentage of ' \
    'ambigious amino acids allowed in a sequence')

    # Define infile
    args = parser.parse_args(sys.argv[1:])
    infile = args.input

    # Set up threshold
    if args.threshold:
        thres = args.threshold
    else:
        thres = '5'

    print('Retrieving data')
    # Retrieve data
    subprocess.run(['python', 'scripts/ncbi_extract.py', '-i', infile, '-o', 'proteins.txt'], check=True)
    print('Data retrieved')

    print(f'Filtering low quality reads with {thres}% threshold')
    # Filter low quality sequences
    subprocess.run(['python', 'scripts/filter.py', '-i', 'proteins.txt', '-o', 
                    'filtered_proteins.txt', '-p', thres])
    print('Reads filtered')

    print('Performing multiple sequence alignment')
    # Perform MSA 
    with open('aligned.txt', 'w') as f:
        subprocess.run(['mafft', '--auto', '--quiet', 'filtered_proteins.txt'], stdout=f, check=True) 
    print('Sequences aligned') 

    print('Running autoencoder')
    # Perform autoencoder
    subprocess.run(['python', 'scripts/autoencoder.py', '-i', 'aligned.txt'])  
    print('AEGIS pipeline complete') 
    

if __name__ == '__main__':
    main()