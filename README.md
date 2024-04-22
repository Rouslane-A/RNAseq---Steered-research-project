# RNAseq---Steered-research-project
## Overview
### Objective
The objective of this repository is to conduct a comprehensive reanalysis of the paper titled 'Single-cell Analysis Reveals Inter- and Intratumour Heterogeneity in Metastatic Breast Cancer' while also performing a novel single-cell analysis using alternative computational methods. The primary goal is to replicate the experimental analyses described in the original paper to assess the reprudicility of the findings. Additionally, this project aims to compare the results obtained from the replication with those derived from the alternative pipelines and software, thereby identifying any additional insights. The repository provides code, data, and documentation to facilitate replication, comparison and interpretation of the analyses conducted. It also includes scripts for the graphical user interface (GUI) which contains the graphs and results from both reanalysis and novel analysis for easier visualization.

### Replication of Original Analysis
The reanalysis of the original paper involves replicating the experimental procedures and computational analyses described in the paper. The code provided in this repository aims to reproduce the results presented in the original publication.
### Novel Analysis
In addition to replicating the original analysis, this project includes a novel single-cell analysis using alternative pipelines and software. the goal of this analysis is to compare the results obtained from different computational approaches and identify any additional insights that may arise from using alternative methodologies.
### Graphical User Interface (GUI)
The graphical user interface developed facilitates the visualization of results from both the reanalysis and the novel analysis. THe GUI provides an interactive platfrom for users to explore and interpret the findings conveniently.

## Contents
**Novel Analysis:** This directory contains all the scripts for the novel analysis, this includes:
- Quality control: Script for assessing the quality of the data
- Filtering: Script for filtering the data based on a defined criteria
- Alignment: Script for aligning the data to the refernce genome
- Count Matrix: Script for genrating the count matrix from aligned data
- Seurat Pipeline: Script detailing the steps involved in using the Seurat package for clustering analysis and visualization.
  
**Reanalysis:** This directory contains all the scripts for the replication of the orginal paper

**Graphical User Interface:** This directory contains scripts for the GUI

## Dependencies
- Access to a high performance computing system particularly ALICE is required.
- The required packages and their corresponding version numbers are listed in the 'requirements.txt' file, ensure these are installed prior to running the scripts.

## How to Use
**Reanalysis**

**Novel Analysis**
- Download the repository
- Navigate to the 'Novel Analysis' directorty
- Navigate to the Quality control directory and run the scripts
- Navigate and run the scripts in the Filtering directory

**Graphical User Interface**

## Why the Project is Useful
This project serves as a valuable resource for:
- Researchers seeking to replicate experimental procedures and analyses from the original study for verification.
- Scientists interested in exploring alternative computational methods in single-cell analysis.
- Professionals looking to identify new insights by comparing results between replicated and novel analyses.
- Individuals who benefit from enhanced accessibility through a graphical user interface for result visualisation.
- Educators and students utilizing it as an educational resource in single-cell analysis.

## Contributors
- BismahGhafoor
- Skycr
- aol3
- Patca26 
- barbaratpferreira
- sca16
- Rouslane-A

## Disclaimer
This project is part of an academic assessment and is designed for educational purposes only.




