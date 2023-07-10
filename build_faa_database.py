import pandas as pd
from Bio import SeqIO

file = open(snakemake.output[0], "w")
for i in range(len(snakemake.input)):
    fasta_sequences = SeqIO.parse(open(snakemake.input[i]),'fasta')
    for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            file.write(">")
            file.write(name)
            file.write("\n")
            file.write(sequence)
            file.write("\n")  
file.close()
