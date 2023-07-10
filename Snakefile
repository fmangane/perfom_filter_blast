import os
SAMPLES = []
rootdir = '/home/manganef/Documents/proteo-phototrophy/snakemake-wd/database/Proteo_sampled'

#rootdir = '/home/manganef/Documents/proteo-phototrophy/snakemake-wd/database/exception'
for rootdir, dirs, files in os.walk(rootdir):
	for subdir in dirs:
		if subdir[0] == "G":
			SAMPLES.append(subdir)

#SAMPLES = ['GCF_905336995.1_ASM90533699v1']
#rule to dezip faa and gff

rule dezipfaa:
	params:
		id = "{id}"
	input:
		faa = "database/Proteo_sampled/{id}/{id}_protein.faa.gz"

	output:
		faa = "database/Proteo_sampled/{id}/{id}_protein.faa"

	shell:
		"gzip -f -d -k {input.faa} > {output.faa}"
rule dezipgff:
	params:
		id = "{id}"
	input:
		gff = "database/Proteo_sampled/{id}/{id}_genomic.gff.gz"

	output:
		gff = "database/Proteo_sampled/{id}/{id}_genomic.gff"

	shell:
		"gzip -f -d -k {input.gff} > {output.gff}"


#rule to create a database for a specie
rule makedb:
	params:
		id = "{id}"
	input:
		"database/Proteo_sampled/{id}/{id}_protein.faa"
	output:
		"database/Proteo_sampled/{id}/{id}_protein.faa.pin"
	shell:
		"makeblastdb -in {input} -dbtype prot > {output}"

#run blast against the new database
rule blastp :
	params:
		id = "{id}"
	input:
		index = "database/Proteo_sampled/{id}/{id}_protein.faa.pin",
		faa = "database/Proteo_sampled/{id}/{id}_protein.faa",
		dshi = "PGC_DShi.fsa"
	output:
		"blasts/blast{id}.txt"
	shell:
		"blastp -db {input.faa} -query {input.dshi} -outfmt \"6 std qcovs\" > {output}"


#rule to create a file with the best_hits using blast_output_manager python script
rule manage_blast :
	params:
		id = "{id}"
	input:
		"blasts/blast{id}.txt",
		"table2.csv",
		"database/Proteo_sampled/{id}/{id}_genomic.gff"
	output:
		"best_hits/best_hits{id}.csv"
	script:
		"scripts/blast_output_manager.py"

#rule to put all the PGC corresponding sequences of a same specie in a file
rule PGC_fasta_collector :
	params:
		id ="{id}"
	input:
		"best_hits/best_hits{id}.csv",
		"database/Proteo_sampled/{id}/{id}_protein.faa"
	output:
		"pgcfolder/PGC{id}.fsa"
	script:
		"scripts/SseqPGC_fasta_collector.py"


#rule to put all the available hit files in the same file

rule get_all_best_hits :
	input:
		csv = expand("best_hits/best_hits{sample}.csv", sample=SAMPLES),
		fsa = expand("pgcfolder/PGC{sample}.fsa", sample=SAMPLES)
	output:
		"all_best_hits.csv"
	shell:
		"awk '(NR == 1) || (FNR > 1)' {input.csv} > {output}"


#rule to put the faa files in the same file
rule get_all_faa_together :
    input:
        expand("database/Proteo_sampled/{sample}/{sample}.faa", sample=SAMPLES)
    output:
        "all_faa.fsa"
    script:
        "scripts/build_faa_database.py"

#rule to get all the hit sequences of a gene in a file
rule get_all_hitseqs :
	params:
		id = "{id}"
	input:
		"all_best_hits.csv",
		"all_faa.fsa"
	output:
		"hitseqfolder/hitseq_{id}.fsa"
	script:
		"scripts/genes_fasta_collector.py"

#rule to create a table  with the count of PGC genes
#rule get_table :
	#input: 
		#"database/list_with_taxonomy.csv",
		#"table2.csv",
		#"all_best_hits.csv"
	#output:
		#"summary_table.csv"
	#script:
		#"scripts/summary_table_maker.py"	
