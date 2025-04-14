# RecombinationHIV

This repository contains the scripts and files for the analysis presented in 'Estimates of HIV-1 within-host recombination rates across the whole genome' Longley et al. The analysis presented in this study considers the within-host measured recombination rate of HIV with sequencing data derived from samples of hundreds of untreated infections in Sub-Saharan Africa. To infer recombination rates, the RATS-LD method from Romero and Feder 2024. The analysis is split into two sections: first, a simulation approach. This looks at sequences generated from SLiM, from 30 simulated infections with genomes 1000bps in length, and 5 distinct 'true' recombination rates are tested on each set of simulated infections, resulting in a total of 150 alignments. Second, we apply the method to the SSA sequencing data collected by the University of Washington as a part of 3 studies of serodifferent transmission pairs, and sequenced by PANGEA. For ethical reasons, it is not possible to publicy share the raw sequencing data. Raw reads and additional metadata can be accessed by application to the PANGEA consortium (www.pangea-hiv.org.).

### Simulations
Pipeline for simulation analysis: 
#### Simulated sequence data
This folder contains 150 fasta files containing the 

Scripts:
1. Slim file
2. Convert .txt to .fasta
3. Identify diverse sites
4. Measure linkage disequilibrium 
5. Estimate recombination 

### Partners 
Prior to the scripts, the raw bam reads must be processsed with Phyloscanner to generate fasta files for each window. 
1. Linkage estimation
2. Filter dataset
3. Produce results and figures
