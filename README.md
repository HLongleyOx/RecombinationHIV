# HIV-1 Within-Host Recombination Rate Analysis

This repository contains the analysis pipeline and supporting files for our study "Estimates of HIV-1 Within-Host Recombination Rates Across the Whole Genome" (Longley et al.). Our research provides comprehensive recombination rate estimates using sequencing data from hundreds of untreated HIV-1 infections in Sub-Saharan Africa.
The analysis employs the RATS-LD method (Romero and Feder 2024) and is structured in two main components:

Simulation Studies: We evaluated method accuracy using SLiM-generated sequences from 30 simulated infections with 1000bp genomes, testing five distinct "true" recombination rates across 150 total alignments.
Empirical Analysis: We applied our validated approach to HIV sequencing data collected by the University of Washington from three studies of serodifferent transmission pairs, sequenced through the PANGEA consortium.

Due to ethical considerations, the raw sequencing data cannot be shared publicly. Researchers can access raw reads and additional metadata through application to the PANGEA consortium (www.pangea-hiv.org).

### Simulations
Pipeline for simulation analysis: 
#### Simulated sequence data
This folder contains 150 fasta files containing the output from SLiM which has been converted from .txt to .fasta. 

#### Scripts:
1. Slim file 
2. Convert .txt to .fasta
   output: Fasta files containing alignments found in the simulated sequence data folder 
4. Identify diverse sites
   Input:
   Output:
6. Measure linkage disequilibrium
   Input:
   Output:
8. Estimate recombination
   Input:
   Output: 

### Partners 
Prior to the scripts, the raw bam reads must be processsed with Phyloscanner to generate fasta files for each window. Phyloscanner default settings are applied, and the set of references in included here as a datafile. Windows of length 750bps are applied to PacBio data and 250bps for Illumina data. Windows over lap by 740bps (or 240bps in case of Illumina). Due to the size of the datasets it was necessary to utilise batch jobs on a HPC environment. 

1. Linkage estimation
   
3. Filter dataset
4. Produce results and figures
