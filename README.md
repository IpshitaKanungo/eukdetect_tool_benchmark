# Benchmarking Tools for Detecting Eukaryotic Sequences in Metagenomic Data
Metagenomic datasets from host-associated environments often contain unwanted host DNA that can distort microbial analyses. Detecting and filtering these eukaryotic (host) sequences is therefore an essential preprocessing step.
This project provides a **benchmarking pipeline** to evaluate how well existing tools identify eukaryotic sequences in metagenomic assemblies—especially when a reference host genome is unavailable.
## Overview
- **Pipeline (Snakemake workflow)**  
Located in the root directory — includes:  
  - `Snakefile`- workflow  
  - `config.yaml`- dataset paths and tool parameters  
  - `environment.yaml`- conda environment  
  - `envs/` - tool-specific environments (eukrep.yaml, tiara.yaml, metawrap.yaml)  
  - `rules/`- Snakemake rules for each step
- **Post-processing scripts:** `data_post-processing`  
    Contains analysis scripts and Jupyter notebooks (in HTML format) for results of benchmarking, computation and visualization
- **Euk & Non-euk FASTA segregation scripts:** `euk_non-euk_segregation_scripts`  
    Contains Python scripts (tiara_split.py, kraken2_split.py) for separating eukaryotic and non-eukaryotic contigs after classification.
## Tools Evaluated
- **Kraken2**
- **EukRep**
- **Tiara**
- *(Baseline)* **Bowtie2** – used when a high-quality host reference genome is available.
## Methodology
1. The **CAMI II Toy Mouse Gut** dataset is used as the test case.
2. Host contamination is simulated by adding synthetic reads from the mouse genome.
3. Reads are assembled into contigs, classified with each tool, and subjected to binning and refinement.
4. Performance is evaluated against ground truth host mappings in terms of:
   - Effect on assembly quality
   - Accuracy in distinguishing eukaryotic vs. non-eukaryotic contigs
   - Influence on final bin quality
## Purpose
This pipeline helps benchmark and compare tools for detecting eukaryotic contamination in metagenomic datasets—both with and without a known host genome
