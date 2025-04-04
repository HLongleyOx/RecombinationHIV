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
    
        for i in range(n):
            seqs = []
            for j in range(m):
                seqs.append(sequences[j][i])

            no_gap_seq = [s for s in seqs if s!="-"]
            a_count = no_gap_seq.count("A")
            c_count = no_gap_seq.count("C")
            g_count = no_gap_seq.count("G")
            t_count = no_gap_seq.count("T")
            count_all = (a_count)*(a_count-1)+(c_count)*(c_count-1)+(g_count)*(g_count-1)+(t_count)*(t_count-1)

            n_len = len(no_gap_seq)*(len(no_gap_seq)-1)
            if len(no_gap_seq)>10:
                div = (n_len-count_all)/n_len
                diversity_df_temp=pd.DataFrame([div, i, date], index=["diversity","site","date"]).T
                diversity_df = pd.concat([diversity_df,diversity_df_temp])
                


        
            if count_occurrences_items(seqs, "-")*(1/m)<0.05 and m>10:
                          
                counts=count_occurrences(seqs)
                counts_all = list(counts.values())
                counts_bases = list(counts.keys())
                if len(counts_all)>1:
                    MAF = counts_all[1]*(1/m)
                    nucleotide = counts_bases[1]
                    if MAF>0.05:
                    
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
   

            
        df=pd.DataFrame([date, siteA, siteB, ddash,baseA,baseB,pA, pB,N, ld], index=["date","siteA","siteB","ddash","baseA","baseB" , "pA","pB", "N", "ld"]).T    
        return(df)
    
    else: 
        return([])

 

        
if __name__ == "__main__":
    couple = sys.argv[1]
    file=sys.argv[2]
    window_len = sys.argv[3]
    print(file)
    window = (file.split("Window_"))[1].split(".fasta")[0]
    window_start = int(window.split("_")[0])
    window_end = int(window.split("_")[2])
    filewindow = "AlignedReadsInWindow_" + str(window_start) + "_to_" + str(window_end) +".fasta"
    
    fastafile = "/well/fraser/users/hqh588/pacbio_updated/data_recipient/data/"+str(couple)+"/phyloscanner_output_750_window/AlignedReads/" + filewindow
    sample_sequences = []

    with open(fastafile) as handle:
        for index, record in SeqIO.FastaIO.SimpleFastaParser(handle):
            if str(couple) in index:
                count = int(index.split("_")[6])
                for x in range(count):
                    #add sequence the number of times it appears
                    ind=index.replace('count_'+str(count), "count_"+str(x))
                    sample_sequences.append([ind,record])

            elif "HXB2" in index:
                hxb2_sequence=record
                
 #find positions of gaps: 
    hxb2_gaps_pos = [idx for idx, item in enumerate(hxb2_sequence) if '-' in item]                         
    sequences_nogaps = []
    for sequence in sample_sequences:
        seq = sequence[1]
        ind = sequence[0]
        sequences_gapfree = ("".join([char for idx, char in enumerate(seq) if idx not in hxb2_gaps_pos]))
        sequences_nogaps.append([ind, sequences_gapfree])



    dates = list(set([a.split("_")[2] for a,b in sequences_nogaps if "HXB2" not in a]))
                                   
    sites_diverse, diversity = find_sites(sequences_nogaps, dates)
    sites_diverse['date_str']=sites_diverse['date']
    sites_diverse['date'] = pd.to_datetime(sites_diverse['date'])
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
        dt = df["date_str"]
        df_othersites = df_othersites[df_othersites["date_str"]==dt]
    
        m = df_othersites.shape[0]
    
        for j in range(m):
            dfB = df_othersites.iloc[j]
            baseB = dfB["first.nucleotide"]
            siteB = dfB["site"]


            linkage_output=linkage(sequences_nogaps, dt, siteA, siteB, baseA,baseB)
    
            if len(linkage_output)>0:
                df_linkage = pd.concat([df_linkage, linkage_output])

        
    df_linkage["couple"]=couple
    df_linkage["window"]=window
    df_linkage["type"]="recip"
    
    diversity["window"]=window
    diversity["couple"]=couple
    diversity["type"]="recip"
    

    if not os.path.exists("/well/fraser/users/hqh588/pacbio_updated/output/diversity/diversity_recipient_" + str(window_len) + "_"+couple+".csv"):
        diversity.to_csv("/well/fraser/users/hqh588/pacbio_updated/output/diversity/diversity_recipient_" + str(window_len) + "_"+couple+".csv",index=False, header=False)
    else:
        diversity.to_csv("/well/fraser/users/hqh588/pacbio_updated/output/diversity/diversity_recipient_" + str(window_len) + "_"+couple+".csv", mode='a', index=False, header=False)


    if not os.path.exists("/well/fraser/users/hqh588/pacbio_updated/output/linkage_recombination/linkage_recipient_" + str(window_len) + "_"+couple+".csv"):
        df_linkage.to_csv("/well/fraser/users/hqh588/pacbio_updated/output/linkage_recombination/linkage_recipient_" + str(window_len) + "_"+couple+".csv",index=False, header=False)

    else:
        df_linkage.to_csv("/well/fraser/users/hqh588/pacbio_updated/output/linkage_recombination/linkage_recipient_" + str(window_len) + "_"+couple+".csv", mode='a', index=False, header=False)



