import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys

#Segregating sites
def parse_segregating_file(output_file):
    mutations_data = []
    genomes_data = []

    with open(output_file, 'r') as file:
        lines = file.readlines()

        # Parse mutations data
        mutations_section_index = lines.index('Mutations:\n') + 1
        genomes_section_index = lines.index('Genomes:\n')
        
        for line in lines[mutations_section_index:genomes_section_index]:
            fields = line.strip().split()
            mutation_id = int(fields[0])
            nucleotide = fields[9]
            position = int(fields[3])  
            mutations_data.append((mutation_id, nucleotide, position))

        # Parse genomes data
        for line in lines[genomes_section_index + 1:]:        
            fields = line.strip().split()
            genome_id = int(fields[0].split(':')[1])
            sequence = ''.join(fields[1:])
            genomes_data.append((genome_id, fields[2:]))
    
    return mutations_data, genomes_data
            
  
#Non segregating sites
def parse_fixed_file(output_file):
    mutations_data = []

    with open(output_file, 'r') as file:
        lines = file.readlines()

        # Parse mutations data
        mutations_section_index = lines.index('Mutations:\n') + 1
        
        for line in lines[mutations_section_index:]:
            fields = line.strip().split()
            mutation_id = int(fields[0])
            position = int(fields[3])
            nucleotide = fields[9]
            mutations_data.append((nucleotide, position))
            
    return mutations_data


def read_fasta(file_path):
    sequence = list()
    with open(file_path) as handle:
        for index, record in SeqIO.FastaIO.SimpleFastaParser(handle):
            sequence.append(record)

    return sequence

def write_sequences_to_fasta(sequences, output_file, time):
    with open(output_file, "w") as f:
        for idx, seq in enumerate(sequences):
            seq_obj = Seq(seq[0])
            seq_record = SeqRecord(seq_obj, id=f"sequence_{time}_{idx}", description="")
            SeqIO.write(seq_record, f, "fasta")

# Example usage
file_path = "/well/fraser/users/hqh588/pacbio_updated/simulations_slim/HXB2.fasta"
sequence = read_fasta(file_path)[0]

ancestral_sequence = []
#make the sequence into a list format 
for i in range(len(sequence)):
    ancestral_sequence.append([sequence[i], i])
    
    
     
output_file = 'gen_' + str(sys.argv[1]) + '_seg.txt'  # Specify the path to your SLiM output file
output_directory = '/well/fraser/users/hqh588/pacbio_updated/simulations_slim/simulations/simulation_rep' + str(sys.argv[2]) + "_rho" + str(sys.argv[3])   # Specify the directory to save the FASTA files
os.chdir(output_directory)
mutations_data, genomes_data = parse_segregating_file(output_file)
fixed_data = parse_fixed_file('gen_' + str(sys.argv[1])+ '_fixed.txt')
 
    
#function to identify positions of segregating sites
n = len(sequence)
m = len(genomes_data)
genomes_all = list()

#If position is not in segrating or fixed, then it is the same as the original 


for j in range(m):
    segregating_muts = genomes_data[j][1] 
    genome = ""                

    muts_numeric = [int(a) for a in segregating_muts]
    segregating_allele = [[b,c] for a,b,c in mutations_data if a in muts_numeric]    
    
    for i in range(len(ancestral_sequence)):
    #segregating first
        seg = [a for a,b in segregating_allele if i==b]
        fix = [a for a,b in fixed_data if i==b]
        anc = [a for a,b in ancestral_sequence if i==b]
        if len(seg)>0:
        #Add onto sequence 
            genome = genome + seg[0]
        elif len(fix)>0:
            genome = genome + fix[0]
        else: 
            genome = genome + anc[0]
                
    genomes_all.append([genome])


# Example usage
time = int(sys.argv[1]) - 50000 
output_file = "/well/fraser/users/hqh588/pacbio_updated/simulations_slim/sequences/output_" + str(sys.argv[2]) + "_rho" + str(sys.argv[3])+ "_" + str(sys.argv[1])  + ".fasta"
write_sequences_to_fasta(genomes_all, output_file, time)
