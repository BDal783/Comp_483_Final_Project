# **COMP483 Project: Machine Learning for Predicting Viral Protein Mutational Landscapes**


## Highlights
  * Users input a template text file to pull FASTA files for viral proteins from NCBI
  * Pipeline creates a multiple sequence alignment using MAFFT
  * Utilizes an autoencoder for unsupervised machine learning and predicts mutational hotspots in viral protein sequences
  * Generates tables and figures to summarize predicted hotspots for analysis

## Overview

In this project, we adapt the AEGIS (AutoEncoder-driven Genomic Insight System) framework, developed by the 
[Miller Lab](https://github.com/mccainwa/Dr.Miller-Lab-MP/blob/main/Instructions.md), 
to streamline the prediction of mutational hotspots in viral proteins. 

Our goal is to build a command-line tool that retrieves FASTA sequences for a 
specified viral protein, performs multiple sequence alignment, and computes 
mutation probabilities at each amino acid position. The tool outputs positions 
most likely to mutate, along with the most probable amino acid substitutions.

This pipeline is designed to simplify large-scale viral sequence analysis and 
provide a fast, reproducible method for identifying regions of evolutionary 
variability. For more information, see the [project wiki](https://github.com/BDal783/Comp_483_Final_Project/wiki)

## Methods
The main script, aegis.py, runs the full pipeline by calling other scripts via subprocess. The general steps the framework implements are as follows:

1. Retrieving viral protein sequences from NCBI: Using the parameter specified in proteinSearch.txt, the pipeline uses Bio.Entrez to search NCBI for the desired viral protein sequences and save them to the FASTA proteins.txt.
2. Filter low-quality reads: Some sequences have a high proportion of ambiguous amino acids (X's), indicating a poor quality read. The pipeline will remove reads that surpass the threshold for the percentage of ambiguous amino acids in the sequence (Ex: threshold = 5, sequences with 5% or greater ambiguous amino acids are removed).
3. Perform MSA using MAFFT: The filtered sequences are aligned using MAFFT. MAFFT will automatically adjust its methods depending on the number of sequences to align to favor speed or accuracy.
4. Format sequences and run the autoencoder: The aligned sequences, stored in the FASTA aligned.txt, are converted into a dataframe format to be fed into the machine learning pipeline. One-hot encoding is then applied to the sequences to convert the amino acids into numerical vectors. The autoencoder model is trained using the TensorFlow package; Monte Carlo dropout inference and Markov transition modeling are performed to factor in model uncertainty and calculate transition probabilities, respectively.
5. Tables and figures: Tables summarizing the top predicted mutation hotspots (local and global positions, amino acid mutations, mutation probabilities) are saved as CSV files. Scatterplots, bar charts, and heatmaps summarizing mutation hotspot probability, Markov transitions, and entropy and variance are saved to the /figures folder.

## Results


## Usage Instructions 
### Setup
This project utilizes a conda virtual environment to ensure the AEGIS pipeline runs smoothly with the correct package versions. Follow the steps below:

Clone repo:
```
git clone https://github.com/BDal783/Comp_483_Final_Project.git
cd Comp_483_Final_Project
```

Install and activate miniconda (Only once):
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/miniconda3/bin/activate
```

Create environment (only once, unless conda_env.yaml changes):
```
conda env create -f conda_env.yaml
```
### Running the Pipeline
Activate virtual environment (every session):
```
conda activate aegis
```
To run the pipeline, edit the template text file called "proteinSearch.txt" with the search parameters for your desired viral protein. You will enter your email address, protein search term, the number of sequences to extract, and the start and end publication dates (format YYYY/MM/DD) to extract sequences from NCBI. When deciding the number of sequences, keep in mind that the pipeline will drop low-quality sequences with a high percentage of ambiguous amino acids, based on the quality threshold defined in aegis.py.

Once you have filled out the template text file, you can run the full pipeline by running this command:
```
python aegis.py -i proteinSearch.txt
```

Deactivate virtual environment (optional)
```
conda deactivate
```

## Authors
This project was worked on by [Brendon Dal](https://github.com/BDal783), [Leah Briscoe](https://github.com/leahbriscoe830), and [Jimmy Capecci](https://github.com/jcapecci09). We are master's students studying bioinformatics at Loyola University Chicago, and we have been tasked with improving the AEGIS framework. This project was given to us by the [Miller Lab](https://wmiller6.sites.luc.edu/assets/GroupMembers.html) to be completed as our final project in COMP483 Computational Biology.

## Special Thanks
We would like to express our appreciation to the members of the Miller lab for providing us with this opportunity and guiding us throughout this project. We would also like to thank our professor, Dr. Heather Wheeler, and our COMP483 classmates for their support and feedback on this project.
