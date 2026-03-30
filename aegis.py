

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

    # Retrieve data
    subprocess.run(['python', 'ncbi_extract.py', '-i', infile, '-o', 'proteins.txt'], check=True)

    # Filter low quality sequences
    subprocess.run(['python', 'filter.py', '-i', 'proteins.txt', '-o', 'filtered_proteins.txt'])

    # Start time of MAFFT
    start_time = time.perf_counter()

    # Perform MSA 
    with open('aligned.txt', 'w') as f:
        subprocess.run(['mafft', '--auto', 'filtered_proteins.txt'], stdout=f, check=True)      

    # find time MAFFT took to run
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print("MAFFT Time: " + str(elapsed_time))

if __name__ == '__main__':
    main()