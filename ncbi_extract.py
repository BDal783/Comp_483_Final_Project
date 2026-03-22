"""Pulls fasta sequences from NCBi based on input file

"""

from Bio import Entrez
from Bio import SeqIO


# Open example file and read in terms
with open('example_input.txt', 'r') as f:
    email = f.readline()
    search_term = f.readline()
    seq_count = f.readline()
    start_date = f.readline()
    end_date = f.readline()

# Set email
Entrez.email = email

# Search for proteins
searchResultHandle = Entrez.esearch(db = "protein", term = search_term, 
                                    retmax = seq_count, idtype = "protein", 
                                    datetype = "pdat", mindate = start_date, 
                                    maxdate = end_date)

# Find id lists
searchResult = Entrez.read(searchResultHandle)           
ids = searchResult["IdList"]                            

# Find fastas based on ids
handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text") 
records = SeqIO.parse(handle, "fasta")

# deduplicates sequences
seen = set()
deduplicated = []

for record in records:
    seq = str(record.seq)
    if seq not in seen:
        deduplicated.append(record)
        seen.add(seq)
