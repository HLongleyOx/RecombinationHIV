#!/bin/python


""" author: Harriet Longley
    purpose: find diverse sites and calculate linkage disequilibrium for recombination rate""" 

from Bio import SeqIO
import os
import pandas as pd
import csv
import sys
from collections import Counter
from collections import defaultdict
import random

#function to identify diverse positions at each 

def count_occurrences(input_list):
    counts = defaultdict(int)  # Automatically initializes each key with 0

    for item in input_list:
        counts[item] += 1

    return counts


def count_occurrences_items(input_list, item):
    return input_list.count(item)


def find_sites(sample_sequences, dates):
        
    diversity_df = pd.DataFrame()
  
    df_return = pd.DataFrame(columns=["site","date","nucleotide"])
    for date in dates:

        sequences = [b for a, b in sample_sequences if date in a]
        
        n = len(sequences[0])
        m = len(sequences)
        print(n, m)
    
        for i in range(n):
            seqs = []
            for j in range(m):
                seqs.append(sequences[j][i])

            a_count = seqs.count("A")
            c_count = seqs.count("C")
            g_count = seqs.count("G")
            t_count = seqs.count("T")
            count_all = (a_count)*(a_count-1)+(c_count)*(c_count-1)+(g_count)*(g_count-1)+(t_count)*(t_count-1)
 
            n_len = len(seqs)*(len(seqs)-1) 
            print(len(seqs))
            div = (n_len-count_all)/n_len
            diversity_df_temp=pd.DataFrame([div, i, date], index=["diversity","site","date"]).T
            diversity_df = pd.concat([diversity_df,diversity_df_temp])

                          
            counts=count_occurrences(seqs)
            counts_all = list(counts.values())
            counts_bases = list(counts.keys())
            if len(counts_all)>1:
                MAF = counts_all[1]*(1/m)
                nucleotide = counts_bases[0]
                if MAF>0.05 and len(counts_all)<3:
                    df_temp = pd.DataFrame([i, date, nucleotide], index=["site","date","nucleotide"]).T
                    df_return=pd.concat([df_return, df_temp])

    
    #Filter so that there are at least 2 mutations 
    grouped_site = df_return.groupby("site")

    # Use the 'filter()' method to filter groups with at least 'min_observations' observations
    filtered_df = grouped_site.filter(lambda x: len(x) >= 2)

    return(filtered_df, diversity_df)

        
    
    
def linkage(sequences_all, date, siteA, siteB, baseA, baseB):
    
    #Filter sites
    
    sequences = [b for a, b in sequences_all if date in a]
    
    seqA = []
    seqB = []
    seqC = []
    for seqs in sequences:
        seqA.append(seqs[siteA])
        seqB.append(seqs[siteB])
        seqC.append(seqs[siteA]+seqs[siteB])
     
    if seqA.count(baseA)==0 or seqB.count(baseB)==0:
        return([])
        
    N = len(seqA)

    pA = seqA.count(baseA)/len(seqA)
    pB = seqB.count(baseB)/len(seqB)
    pa = 1-pA
    pb = 1-pB
    pAB = seqC.count(baseA + baseB)/len(seqC)
    ld = pAB - pA*pB
    for i in range(len(seqC)):
        read=seqC[i]
        if read[0]!=baseA and read[1]!=baseB:
            seqC[i] = "NN"
        elif read[0]==baseA and read[1]!=baseB:
            seqC[i] = baseA + "N"
        elif read[0]!=baseA and read[1]==baseB:
            seqC[i] = "N" + baseB
    
    
    haploFreq=pAB
    counts = len(seqC)
    A = seqA.count(baseA)
    a = len(seqA)-A
    B = seqB.count(baseB)
    b = len(seqB)-B
        
    basea="N"
    baseb="N"
    
    AB = float(seqC.count(baseA + baseB))
    aB = float(seqC.count(basea + baseB))
    Ab = float(seqC.count(baseA + baseb))
    ab = float(seqC.count(basea + baseb))
    

    
    if AB>0 and aB>0 and Ab>0 and ab>0:
        maxF = max(-1*(pA*pB), -1*(pa*pb))
        minF = min((pA*pb),(pa*pB))

        if ld<0:
            ddash = ld/maxF
        elif ld>0:
            ddash = ld/minF
   
        else: 
            ddash=0
            
        df=pd.DataFrame([date, siteA, siteB, ddash, baseA,baseB, pA, pB,pAB, N,ld], index=["date","siteA","siteB","ddash","baseA","baseB","pA","pB","pAB", "N","ld"]).T    
        return(df)
    
    else: 
        return([])

 

        
if __name__ == "__main__":
    os.chdir("/well/fraser/users/hqh588/pacbio_updated/simulations_slim")
    sim=sys.argv[1]
    fastafile = "sequences/sequences_rho" + str(sys.argv[2])+ "_" + str(sim) +  ".fasta"
    sample_sequences = []
    
    with open(fastafile) as handle:
        for index, record in SeqIO.FastaIO.SimpleFastaParser(handle):
            sample_sequences.append([index,record])
    
    dates=["_" + str(x) + "_" for x in range(0, 601,50)]
    sites_diverse, diversity = find_sites(sample_sequences, dates)
    
    sites_diverse['date_str']=sites_diverse['date']
    sites_diverse['date'] = sites_diverse['date']
    sites_diverse = sites_diverse.sort_values(by=["site","date"])
    mask = sites_diverse['date'] == sites_diverse.groupby('site')['date'].transform('min')
    sites_diverse_first=sites_diverse[mask].rename(columns={'nucleotide': 'first.nucleotide', 'date': 'first.date'}).drop("date_str", axis=1)
    sites_diverse = pd.merge(sites_diverse, sites_diverse_first, on="site")  
    
    
    n = sites_diverse.shape[0]
    df_linkage = pd.DataFrame()
    
    for i in range(n):
        df = sites_diverse.iloc[i]
        baseA = df["first.nucleotide"]
        siteA = df["site"]
        df_othersites = sites_diverse.iloc[(i+1):]
        dt = df["date"]
        df_othersites = df_othersites[df_othersites["date"]==dt]
    
        m = df_othersites.shape[0]
    
        for j in range(m):
            dfB = df_othersites.iloc[j]
            baseB = dfB["first.nucleotide"]
            siteB = dfB["site"]


            linkage_output=linkage(sample_sequences, dt, siteA, siteB, baseA,baseB)
    
            if len(linkage_output)>0:
                df_linkage = pd.concat([df_linkage, linkage_output])
    
    df_linkage["ID"] = sim
    df_linkage["rho"] = sys.argv[2]
    df_linkage.to_csv("linkage/linkage_" +str(sim)+ "_rho" + str(sys.argv[2]) + ".csv")

