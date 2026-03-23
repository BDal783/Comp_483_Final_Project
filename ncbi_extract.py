"""Pulls fasta sequences from NCBI based on input file

Input file must have the following format
Email
NCBI search term
Number of sequences
min date
max date
"""

from Bio import Entrez
from Bio import SeqIO
import argparse
import sys


def main():
    # Set up parser for command line interface
    parser = argparse.ArgumentParser(description='Pulls fasta sequences from NCBI')
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)

    # Define infile and outfile
    args = parser.parse_args(sys.argv[1:])
    infile = args.input
    outfile = args.output

    # Open example file and read in terms
    with open(infile, 'r') as f:
        email = f.readline().strip()
        search_term = f.readline().strip()
        seq_count = int(f.readline().strip())
        start_date = f.readline().strip()
        end_date = f.readline().strip()

    # Set email
    Entrez.email = email

    # Search for proteins
    searchResultHandle = Entrez.esearch(db = "protein", term = search_term, 
                                        retmax = seq_count, datetype = "pdat", 
                                        mindate = start_date, maxdate = end_date)

    # Find id lists
    searchResult = Entrez.read(searchResultHandle)           
    ids = searchResult["IdList"]                         

    # Find fastas based on ids
    handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text") 
    records = SeqIO.parse(handle, "fasta")

    # intialize  set and list to deduplicate sequences
    seen = set()
    deduplicated = []

    # remove any duplicates from sequences
    for record in records:
        seq = str(record.seq)

        #  If seq  is not in the set add it to the deduplicated list
        # and the set
        if seq not in seen:
            deduplicated.append(record)
            seen.add(seq)
    
    # Write fastas to mulit-fasta file
    with open(outfile, 'w') as f:
        SeqIO.write(deduplicated, f, "fasta")


if __name__  == '__main__':
    main()
