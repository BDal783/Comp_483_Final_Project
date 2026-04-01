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
variability. For more information, see the [project wiki](https://github.com/BDal783/Comp_483_Final_Project/wiki)

## Results



## Usage instructions 
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

Create environment (only once, unless conda_env.yaml changes)
```
conda env create -f conda_env.yaml
```

Activate virtual environment (every session):
```
conda activate imlabtools2 
```
Run the pipeline with example data:
```
python aegis.py -i proteinSearch.txt
python aegis.py -i proteinSearch.txt
```

Deactivate virtual environment (optional)
```
conda deactivate
```

## Authors
This project was worked on by and [Brendon Dal](https://github.com/BDal783), [Leah Briscoe](https://github.com/leahbriscoe830), and [Jimmy Capecci](https://github.com/jcapecci09). We are masters students studying bioinformatics at Loyola University Chicago, who have been tasked with improving the AEGIS framework. This project was given to us by the [Miller Lab](https://wmiller6.sites.luc.edu/assets/GroupMembers.html) to be compleetd as our final project in Comp483.
