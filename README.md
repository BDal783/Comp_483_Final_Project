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

1. Retrieving viral protein sequences from NCBI
2. Filter low-quality reads
3. Perform MSA using MAFFT
4. one-hot encoding (incude gaps and ambigious amino acids into encoding scheme)
5. run autoencoder and include Monte Carlo dropout inference and Markov transition modeling 
7. Output figures

## Outputs

* heatmap_probs_top_hotspots.png: Shows mutational hotspots and probabiliy of mutating to each amino acid
* entropy_strip.png: High entropy peaks represent areas where the Autoencoder’s "accuracy" is lowest. This tells you which parts of the protein the model finds complex or poorly represented in the training data.
* 


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
To run the pipeline, edit the template text file called "proteinSearch.txt" with the search parameters for your desired viral protein. You will enter your email address, protein search term, the number of sequences to extract, and the start and end publication dates (format YYYY/MM/DD) to extract sequences from NCBI. When deciding the number of sequences, keep in mind that the pipeline will drop low-quality sequences with a high percentage of ambiguous amino acids, based on the quality threshold defined in aegis.py. The maximum number of sequences you can extract is 10,000.

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
