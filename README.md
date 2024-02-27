# perfom_filter_blast
A pipeline to perform BLAST on several genomes using one specie as a model then filter the output to keep only the best hits (based on e-value, coverage and bitscore).
The analysis adds annotation informations ( gene name, function, position on genome...) to the hits remaining after filters.
Initial files to be provided : 
_ fasta format of the search protein model 
_ for each genome faa and gff annotations files.

# Usage 
Working directory must be where the Snakefile is.
Folders : database containing each genome and their required files ; scripts with all python scripts in it ; file containing the model ( DShi_PGC.txt in this exemple). 
Then 

'''
snakemake -c 1 
''' 
