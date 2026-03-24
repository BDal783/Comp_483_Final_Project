

import subprocess
import argparse
import sys

def main():

    # Set up parser
    parser = argparse.ArgumentParser(description='Runs the AEGIS pipeline')
    parser.add_argument('-i', '--input', help='input file', required=True)

    # Define infile
    args = parser.parse_args(sys.argv[1:])
    infile = args.input

    # Retrieve data
    subprocess.run(['python', 'ncbi_extract.py', '-i', infile, '-o', 'output.txt'])

if __name__ == '__main__':
    main()