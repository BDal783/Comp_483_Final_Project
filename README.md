# **Comp_483_Final_Project**


## Highlights
  * Successfully identifies mutational hotspots in viral proteins
  * Users input text file prompting software to pull fasta files from NCBI
  * Pipeline performs MSA using [Insert tool]
  * Utilizes an autoencoder in an unsupervised fashion and returns amino acids most likely to mutate

## Overview

In this project we adapt the AEGIS (AutoEncoder-driven Genomic Insight System) framework, developed by the 
[Miller Lab](https://github.com/mccainwa/Dr.Miller-Lab-MP/blob/main/Instructions.md), 
to streamline the identification of mutational hotspots in viral proteins. 

Our goal is to build a command-line tool that retrieves FASTA sequences for a 
specified viral protein, performs multiple sequence alignment, and computes 
mutation probabilities at each amino acid position. The tool outputs positions 
most likely to mutate, along with the most probable amino acid substitutions.

This pipeline is designed to simplify large-scale viral sequence analysis and 
provide a fast, reproducible method for identifying regions of evolutionary 
variability.


## Data Retrieval 
The pipeline utilizes the Entrez module from Biopython to pull FASTA sequences from NCBI. In order, to retrieve sequences the user must input a text file (Ex: proteinSearch.txt)
including the following information: 

Email: {insert email}  
Protein name: {Insert protein}  
Number of sequences to extract: {insert number}  
Start date: {insert date in yyyy/mm/dd format}  
End date: {insert date in yyyy/mm/dd format}  

This process is performed using the ncbi_extract.py script in the following fashion:
```
python ncbi_extract.py -i proteinSearch.txt -o example_fasta.txt
```
However, refer to usage instructions on how to use command-line tool


## Multiple Sequence Alignment (MSA)


## Authors
This project was worked on by and [Brendon Dal](https://github.com/BDal783), [Leah Briscoe](https://github.com/leahbriscoe830), and [Jimmy Capecci](https://github.com/jcapecci09). We are masters students studying bioinformatics at Loyola University Chicago, who have been tasked with improving the AEGIS framework. This project was given to us by the [Miller Lab](https://wmiller6.sites.luc.edu/assets/GroupMembers.html) to be compleetd as our final project in Comp483.


## Usage instructions 
This project utilizes a conda virtual enviornment to ensure the AEGIS pipeline runs smoothly with the correct package versions. The following are instructions on setting up the enviornment and running the AEGIS pipeline.   

Install and activate miniconda: Only needs to be executed once
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/miniconda3/bin/activate
```

Create enviornment: Only needs to be executed once unless contents of conda_env.yaml are altered
```
conda env create -f conda_env.yaml
```

Activate virtual enviornment: Needs to be executed everytime you run pipeline
```
conda activate imlabtools2 
```

Run pipeline with example data:
```
python aegis.py -i proteinSearch.txt
```

Deactivate virtual enviornment:
```
conda deactivate
```