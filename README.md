# **Comp_483_Final_Project**


## Highlights
  * Successfully identifies mutational hotspots in viral proteins
  * Users input text file prompting software to pull fasta files from NCBI
  * Pipeline performs MSA using [Insert tool]
  * Utilizes an autoencoder in an unsupervised fashion and returns amino acids most likely to mutate

## Overview

In this project we adapt the AEGIS framework, developed by the 
[Miller Lab](https://github.com/mccainwa/Dr.Miller-Lab-MP/blob/main/Instructions.md), 
to streamline the identification of mutational hotspots in viral proteins. 

Our goal is to build a command-line tool that retrieves FASTA sequences for a 
specified viral protein, performs multiple sequence alignment, and computes 
mutation probabilities at each amino acid position. The tool outputs positions 
most likely to mutate, along with the most probable amino acid substitutions.

This pipeline is designed to simplify large-scale viral sequence analysis and 
provide a fast, reproducible method for identifying regions of evolutionary 
variability.

## data retrieval 
The pipeline utilizes the Entrez module from Biopython to pull FASTA sequences from NCBI. In order, to retrieve sequences the user must input a text file (Ex: proteinSearch.txt)
including the following information: 

Email: {insert email}
Protein name: {Insert protein}
Number of sequences to extract: {insert number}
Start date: {insert date in yyyy/mm/dd format}
End date: {insert date in yyyy/mm/dd format}

This process is performed using the ncbi_extract_merged.py script in the following fashion:
```
python ncbi_extract_merged.py -i proteinSearch.txt -o muttifasta.txt
```
However, refer to usage instructions on how to use command-line tool


## Authors


## Usage instructions 

Project is work in progress.
[Insert instructions]

