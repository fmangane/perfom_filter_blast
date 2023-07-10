
import pandas as pd
from Bio import SeqIO


accession = pd.read_csv(snakemake.input[0],header=0)
accessionlist = []
for i in range(len(accession)):
    accessionlist.append(accession['sseqid'][i])
file = open(snakemake.output[0], "w")
#parser for proteine.faa
protein = snakemake.input[1]
fasta_sequences = SeqIO.parse(open(protein),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    if name in accessionlist :
        file.write(">")
        file.write(name)
        file.write("\n")
        file.write(sequence)
        file.write("\n")
file.close()


