# **COMP483 Project: Machine Learning for Predicting Viral Protein Mutational Landscapes**


## 💫 Highlights
  * Users input a template text file to pull FASTA files for viral proteins from NCBI
  * Pipeline creates a multiple sequence alignment using MAFFT
  * Utilizes an autoencoder for unsupervised machine learning and predicts mutational hotspots in viral protein sequences
  * Generates tables and figures to summarize predicted hotspots for analysis

## 📑 Overview

In this project, we adapt the AEGIS (AutoEncoder-driven Genomic Insight System) framework, developed by the 
[Miller Lab](https://github.com/mccainwa/Dr.Miller-Lab-MP), 
to streamline the prediction of mutational hotspots in viral proteins. 

Our goal is to build a command-line tool that retrieves FASTA sequences for a 
specified viral protein, performs multiple sequence alignment, and computes 
mutation probabilities at each amino acid position. The tool outputs positions 
most likely to mutate, along with the most probable amino acid substitutions.

This pipeline is designed to simplify large-scale viral sequence analysis and 
provide a fast, reproducible method for identifying regions of evolutionary 
variability. For more information, see the [project wiki](https://github.com/BDal783/Comp_483_Final_Project/wiki)

## ⚙️ Methods
The main script, aegis.py, runs the full pipeline by calling other scripts via subprocess. For more information, visit the [project wiki](https://github.com/BDal783/Comp_483_Final_Project/wiki). The general steps the framework implements are as follows:

1. Retrieving viral protein sequences from NCBI
2. Filter low-quality reads
3. Perform MSA using MAFFT
4. One-hot encoding (include gaps and ambiguous amino acids into the encoding scheme)
5. Run autoencoder and include Monte Carlo dropout inference and Markov transition modeling 
7. Output figures

## ✨ Usage Instructions 
### 🛠️ Setup
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
### ▶️ Running the Pipeline
Activate virtual environment (every session):
```
conda activate aegis
```
To run the pipeline, edit the template text file called "proteinSearch.txt" with the search parameters for your desired viral protein. You will enter your email address, protein search term, the number of sequences to extract, and the start and end publication dates (format YYYY/MM/DD) to extract sequences from NCBI. When deciding the number of sequences, keep in mind that the pipeline will drop low-quality sequences with a high percentage of ambiguous amino acids, based on the quality threshold defined in aegis.py. The maximum number of sequences you can extract is 10,000.

Once you have filled out the template text file, you can run the full pipeline by running this command:
```
python aegis.py -i proteinSearch.txt
```

If you would like to adjust the filtering threshold for the amount of ambiguous amino acids, enter an additional argument with the percent cutoff; the default is set to 5%. As an example, to set the threshold to remove 10% or greater ambiguous sequences, run the following command:
```
python aegis.py -i proteinSearch.txt -p 10
```

Deactivate virtual environment (optional)
```
conda deactivate
```

## 🔍 Outputs

<img width="1595" height="600" alt="Screenshot 2026-05-01 090205" src="https://github.com/user-attachments/assets/7a100a1b-097f-4205-ac01-d1768d741630" />

* heatmap_probs_top_hotspots.png: Shows mutational hotspots and the probability of each amino acid occurring at that position. 

<img width="1595" height="600" alt="Screenshot 2026-05-01 090214" src="https://github.com/user-attachments/assets/758e503d-4ae6-43c1-839b-a7b38021d94e" />

* heatmap_probs_normalized.png: Normalized version of the previous heatmap that sets the top probability to 100%. Highlights the variability at that position. 

<img width="1595" height="600" alt="Screenshot 2026-05-01 141817" src="https://github.com/user-attachments/assets/59f2e67d-3048-4485-8c69-b5a20a66a90e" />

* heatmap_probs_least_variable.png: Shows the sequence positions with the least variability in amino acids. One cell being yellow indicates that only that amino acid shows up in that position across variants.

## 🙏 Authors
This project was worked on by [Brendon Dal](https://github.com/BDal783), [Leah Briscoe](https://github.com/leahbriscoe830), and [Jimmy Capecci](https://github.com/jcapecci09). We are master's students studying bioinformatics at Loyola University Chicago, and we have been tasked with improving the AEGIS framework. This project was given to us by the [Miller Lab](https://wmiller6.sites.luc.edu/assets/GroupMembers.html) as our final project in COMP483 Computational Biology.

## 👏 Special Thanks
We would like to express our appreciation to the members of the Miller lab for providing us with this opportunity and guiding us throughout this project. We would also like to thank our professor, Dr. Heather Wheeler, and our COMP483 classmates for their support and feedback on this project.
