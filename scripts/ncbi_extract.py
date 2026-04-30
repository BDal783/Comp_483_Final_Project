"""
Takes a template txt named proteinSearch.txt as input and extracts the search parameters to 
search the NCBI protein database for viral proteins. Search results saved to proteinResults.txt.

proteinSearch.txt has the following format:

Email:                           Enter email to access Entrez
Protein name:                    Search terms for a viral protein
Number of sequences to extract:  The number of viral protein sequences to retreive from NCBI
Start date:                      The start publication date in the format YYYY/MM/DD
End date:                        The end publication date in the format YYYY/MM/DD

Author: leahbriscoe830 Leah Briscoe and jcapecci09 Jimmy Capecci
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

    # Function to parse proteinSearch.txt for the search parameters
    def protein_search(input: str) -> dict[str, str]:
        """Turns input file into a dictionary based on colon 

        ex: {'Email': 'youremail@gmail.com',.....}

        :param input: A  input file
        :return: Dictionary with contents of file
        """

        # Seperate items from input file and place into a dictionary
        search_terms = {}   
        for line in input:
            param_line = line.split(':')
            search_terms[param_line[0]] = param_line[1].strip()
        return search_terms

    # Read input file and parse search parameters  
    with open(infile, 'r') as i:
        search = protein_search(i)                         # extract search parameters from the input file

    # Protein Sequence Retrieval from NCBI based on search terms
    Entrez.email = search['Email']                         # user email (tell NCBI who you are to access sequences)

    protTerm = search['Protein name']                      # protein name/term
    numSeqs = search['Number of sequences to extract']     # number of seqs to retrieve

    startDate = search['Start date']                       # start date for date range
    endDate = search['End date']                           # end date for date range

    # Create search handle using Entrez.esearch
    searchResultHandle = Entrez.esearch(db = "protein", term = protTerm, 
                                        retmax = numSeqs, idtype = "protein", datetype = "pdat", 
                                        mindate = startDate, maxdate = endDate)                     # entrez search handle with user-specified options
    searchResult = Entrez.read(searchResultHandle)                                       # read in handle with parameters set to user inputs 
    ids = searchResult["IdList"]                                                         # list of IDs created from protein sequences retrieved

    handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text") # protein sequences retrieved by IDs in fasta format
    records = SeqIO.parse(handle, "fasta")                              # record created reading in the handle containing the fasta format protein sequences           


    # Write fastas to mulit-fasta file
    with open(outfile, 'w') as f:
        SeqIO.write(records, f, "fasta")

if __name__  == '__main__':
    main()
    