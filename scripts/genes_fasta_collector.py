import pandas as pd
from Bio import SeqIO
file = open(snakemake.output[0],"w")
df = pd.read_csv(snakemake.input[0],sep =",")
hitlist = []
organismlist = []
for i in range(len(df)):
    if df.loc[i,'qseqid'] == snakemake.params[0] :
        hitlist.append(df.loc[i,'sseqid'])
        organismlist.append(df.loc[i,'Organism'])
dico = {}
protein = snakemake.input[1]
fasta_sequences = SeqIO.parse(open(protein),'fasta')
for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        dico[name]=sequence
for i in range(len(hitlist)):
    if hitlist[i] in  list(dico.keys()):
        file.write(">")
        file.write(hitlist[i])
        file.write(" ")
        file.write(organismlist[i])
        file.write("\n")
        file.write(dico[hitlist[i]])
        file.write("\n")
file.close()
