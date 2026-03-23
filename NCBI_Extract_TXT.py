"""
Updated 3/19/2026 by leahbriscoe830

NCBI_Extract_TXT.py: Takes a template txt named proteinSearch.txt as input and extracts the search parameters to 
search the NCBI protein database for viral proteins. Search results saved to proteinResults.txt.

proteinSearch.txt has the following format:

Email:                           Enter email to access Entrez
Protein name:                    Name of the viral protein
Number of sequences to extract:  The number of viral protein sequences to retreive from NCBI
Start date:                      The start publication date in the format YYYY/MM/DD
End date:                        The end publication date in the format YYYY/MM/DD
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:57:36 2023

@author: blim
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:49:58 2023

@author: oprincewell
"""

# Import necessary packages
import sys                                              # module contains methods and variables for modifying Python's runtime environment
import csv                                              # module implements classes to read and write tabular data in csv format
from Bio import SeqIO                                   # module functioning as an interface to input and output fasta format files
from Bio import Entrez                                  # module to search NCBI for protein sequences with user-specified parameters
from datetime import date                               # module to pull current dates
                                  

# Function to parse proteinSearch.txt for the search parameters
def protein_search(infile):
    search_terms = {}   
    for line in infile:
        param_line = line.split(':')
        search_terms[param_line[0]] = param_line[1].strip()
    return search_terms

# Read input file and parse search parameters
input_file = open('proteinSearch.txt', 'r')           # open input file proteinSearch.txt     
search = protein_search(input_file)                   # extract search parameters from the input file
input_file.close()

# Protein Sequence Retrieval from NCBI based on search terms
Entrez.email = search['Email']                         # user email (tell NCBI who you are to access sequences)

protTerm = search['Protein name']                      # protein name/term
numSeqs = search['Number of sequences to extract']     # # seqs to retrieve

startDate = search['Start date']                       # start date for date range
today = date.today()                                   # current date pulled using datetime module for end date for default date range
endDate = search['End date']                           # end date for date range

# Create search handle using Entrez.esearch
searchResultHandle = Entrez.esearch(db = "protein", term = protTerm, retmax = numSeqs, idtype = "protein", datetype = "pdat", mindate = startDate, maxdate = endDate) # entrez search handle with user-specified options
searchResult = Entrez.read(searchResultHandle)          # read in handle with parameters set to user inputs 
ids = searchResult["IdList"]                            # list of IDs created from protein sequences retrieved

handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text") # protein sequences retrieved by IDs in fasta format
record = handle.read()                                  # record created reading in the handle containing the fasta format protein sequences           

# Write Entrez search results to the output file
output_file = open('proteinResults.txt', 'w')           # open output file proteinResults.txt     
output_file.write(record.rstrip('\n'))                  # each fasta format protein sequence is stripped of the new line character
output_file.close()                                     # close the output file containing the fasta format protein sequences