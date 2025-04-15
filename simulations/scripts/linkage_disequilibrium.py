#!/bin/python

""" 
author: Harriet Longley
purpose: find diverse sites from simulation data and calculate linkage disequilibrium for recombination rate
""" 

from Bio import SeqIO  # For parsing FASTA files
import os
import pandas as pd
import csv
import sys
from collections import Counter
from collections import defaultdict
import random

# Function to count occurrences of each unique item in a list
# Returns a dictionary with items as keys and counts as values
def count_occurrences(input_list):
    counts = defaultdict(int)  # Automatically initializes each key with 0
    for item in input_list:
        counts[item] += 1
    return counts

# Function to count occurrences of a specific item in a list
# This function appears redundant since it just calls the built-in .count() method
def count_occurrences_items(input_list, item):
    return input_list.count(item)

# Function to identify sites with genetic diversity within samples across different timepoints
def find_sites(sample_sequences, dates):
    # Initialize empty dataframes to store results    
    diversity_df = pd.DataFrame()
    df_return = pd.DataFrame(columns=["site", "date", "nucleotide"])
    
    # Process sequences for each sampling date
    for date in dates:
        # Extract sequences for the current date
        sequences = [b for a, b in sample_sequences if date in a]
        
        n = len(sequences[0])  # Number of sites (sequence length)
        m = len(sequences)     # Number of sequences
        print(n, m)
    
        # Analyze each nucleotide position in the sequences
        for i in range(n):
            # Extract the nucleotide at position i from all sequences
            seqs = []
            for j in range(m):
                seqs.append(sequences[j][i])

            # Count occurrences of each nucleotide
            a_count = seqs.count("A")
            c_count = seqs.count("C")
            g_count = seqs.count("G")
            t_count = seqs.count("T")
            
            # Calculate homogeneity sum for diversity calculation
            # Formula counts homozygous pairs (nx * (nx-1)) for each nucleotide x
            count_all = (a_count)*(a_count-1)+(c_count)*(c_count-1)+(g_count)*(g_count-1)+(t_count)*(t_count-1)
 
            # Calculate total possible pairs
            n_len = len(seqs)*(len(seqs)-1) 
            print(len(seqs))
            
            # Calculate nucleotide diversity (proportion of different nucleotides in all pairs)
            div = (n_len-count_all)/n_len
            
            # Store diversity information
            diversity_df_temp = pd.DataFrame([div, i, date], index=["diversity", "site", "date"]).T
            diversity_df = pd.concat([diversity_df, diversity_df_temp])

            # Count occurrences of each nucleotide at this position          
            counts = count_occurrences(seqs)
            counts_all = list(counts.values())
            counts_bases = list(counts.keys())
            
            # If there is variation at this site (more than one nucleotide type)
            if len(counts_all) > 1:
                # Calculate minor allele frequency
                MAF = counts_all[1] * (1/m)
                nucleotide = counts_bases[0]
                
                # Filter for sites with MAF > 0.05 and exactly 2 nucleotide types
                if MAF > 0.05 and len(counts_all) < 3:
                    df_temp = pd.DataFrame([i, date, nucleotide], index=["site", "date", "nucleotide"]).T
                    df_return = pd.concat([df_return, df_temp])

    # Filter so that we only keep sites that appear in at least 2 timepoints
    grouped_site = df_return.groupby("site")
    filtered_df = grouped_site.filter(lambda x: len(x) >= 2)

    return(filtered_df, diversity_df)

# Function to calculate linkage disequilibrium between two sites
def linkage(sequences_all, date, siteA, siteB, baseA, baseB):
    # Filter sequences by date
    sequences = [b for a, b in sequences_all if date in a]
    
    # Extract nucleotides at the two sites
    seqA = []
    seqB = []
    seqC = []  # Combined haplotypes
    for seqs in sequences:
        seqA.append(seqs[siteA])
        seqB.append(seqs[siteB])
        seqC.append(seqs[siteA] + seqs[siteB])
     
    # Return empty list if either of the specified bases is not present
    if seqA.count(baseA) == 0 or seqB.count(baseB) == 0:
        return([])
        
    N = len(seqA)  # Total number of sequences

    # Calculate allele frequencies
    pA = seqA.count(baseA) / len(seqA)  # Frequency of allele A
    pB = seqB.count(baseB) / len(seqB)  # Frequency of allele B
    pa = 1 - pA  # Frequency of non-A alleles
    pb = 1 - pB  # Frequency of non-B alleles
    
    # Calculate haplotype frequency
    pAB = seqC.count(baseA + baseB) / len(seqC)  # Frequency of A-B haplotype
    
    # Calculate raw linkage disequilibrium
    ld = pAB - pA * pB
    
    # Recode haplotypes for counting
    for i in range(len(seqC)):
        read = seqC[i]
        if read[0] != baseA and read[1] != baseB:
            seqC[i] = "NN"  # Neither allele matches
        elif read[0] == baseA and read[1] != baseB:
            seqC[i] = baseA + "N"  # Only first allele matches
        elif read[0] != baseA and read[1] == baseB:
            seqC[i] = "N" + baseB  # Only second allele matches
    
    haploFreq = pAB  # Haplotype frequency
    counts = len(seqC)  # Total count
    
    # Count alleles
    A = seqA.count(baseA)  # Count of allele A
    a = len(seqA) - A      # Count of non-A alleles
    B = seqB.count(baseB)  # Count of allele B
    b = len(seqB) - B      # Count of non-B alleles
        
    basea = "N"  # Placeholder for non-A allele
    baseb = "N"  # Placeholder for non-B allele
    
    # Count haplotypes
    AB = float(seqC.count(baseA + baseB))  # A-B haplotype count
    aB = float(seqC.count(basea + baseB))  # non-A, B haplotype count
    Ab = float(seqC.count(baseA + baseb))  # A, non-B haplotype count
    ab = float(seqC.count(basea + baseb))  # non-A, non-B haplotype count
    
    # Calculate standardized linkage disequilibrium only if all haplotypes exist
    if AB > 0 and aB > 0 and Ab > 0 and ab > 0:
        # Calculate bounds for D' (standardized LD)
        maxF = max(-1 * (pA * pB), -1 * (pa * pb))
        minF = min((pA * pb), (pa * pB))

        # Calculate D' based on sign of LD
        if ld < 0:
            ddash = ld / maxF
        elif ld > 0:
            ddash = ld / minF
        else: 
            ddash = 0
            
        # Create dataframe with results
        df = pd.DataFrame([date, siteA, siteB, ddash, baseA, baseB, pA, pB, pAB, N, ld], 
                          index=["date", "siteA", "siteB", "ddash", "baseA", "baseB", "pA", "pB", "pAB", "N", "ld"]).T    
        return(df)
    
    else: 
        return([])  # Return empty list if any haplotype is missing

# Main execution block
if __name__ == "__main__":
    os.chdir("simulations_slim")  # Change to directory where simulations are saved
    
    # Get simulation parameters from command line arguments
    sim = sys.argv[1]  # Simulation ID
    rho = sys.argv[2]  # Recombination rate parameter
    
    # Construct path to FASTA file
    fastafile = "sequences/sequences_rho" + str(rho) + "_" + str(sim) + ".fasta"
    
    # Load sequences from FASTA file
    sample_sequences = []
    with open(fastafile) as handle:
        for index, record in SeqIO.FastaIO.SimpleFastaParser(handle):
            sample_sequences.append([index, record])
    
    # Define sampling timepoints (dates)
    dates = ["_" + str(x) + "_" for x in range(0, 601, 50)]
    
    # Find diverse sites and calculate diversity metrics
    sites_diverse, diversity = find_sites(sample_sequences, dates)
    
    # Process diverse sites data
    sites_diverse['date_str'] = sites_diverse['date']
    sites_diverse = sites_diverse.sort_values(by=["site", "date"])
    
    # Identify first occurrence of each site and create reference dataframe
    mask = sites_diverse['date'] == sites_diverse.groupby('site')['date'].transform('min')
    sites_diverse_first = sites_diverse[mask].rename(columns={'nucleotide': 'first.nucleotide', 'date': 'first.date'}).drop("date_str", axis=1)
    
    # Merge to add reference information to all occurrences
    sites_diverse = pd.merge(sites_diverse, sites_diverse_first, on="site")  
    
    # Calculate linkage disequilibrium for all pairs of diverse sites
    n = sites_diverse.shape[0]
    df_linkage = pd.DataFrame()
    
    # For each site
    for i in range(n):
        df = sites_diverse.iloc[i]
        baseA = df["first.nucleotide"]  # Reference nucleotide for site A
        siteA = df["site"]              # Position of site A
        dt = df["date"]                 # Date of comparison
        
        # Get all other sites from the same date to compare with
        df_othersites = sites_diverse.iloc[(i+1):]
        df_othersites = df_othersites[df_othersites["date"] == dt]
    
        m = df_othersites.shape[0]
    
        # For each other site
        for j in range(m):
            dfB = df_othersites.iloc[j]
            baseB = dfB["first.nucleotide"]  # Reference nucleotide for site B
            siteB = dfB["site"]              # Position of site B

            # Calculate linkage disequilibrium
            linkage_output = linkage(sample_sequences, dt, siteA, siteB, baseA, baseB)
    
            # Add to results if valid output
            if len(linkage_output) > 0:
                df_linkage = pd.concat([df_linkage, linkage_output])
    
    # Add simulation parameters to results
    df_linkage["ID"] = sim
    df_linkage["rho"] = rho
    
    # Save results to CSV file
    df_linkage.to_csv("linkage/linkage_" + str(sim) + "_rho" + str(rho) + ".csv")
