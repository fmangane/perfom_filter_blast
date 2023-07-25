import os
SAMPLES = []
rootdir = '/database'
for rootdir, dirs, files in os.walk(rootdir):
	if dirs[0] == "G":
		SAMPLES.append(dirs)

#rule to dezip faa and gff

rule dezipfaa:
	params:
		id = "{id}"
	input:
		faa = "database/{id}_protein.faa.gz"

	output:
		faa = "database/{id}_protein.faa"

	shell:
		"gzip -f -d -k {input.faa} > {output.faa}"
rule dezipgff:
	params:
		id = "{id}"
	input:
		gff = "database/{id}_genomic.gff.gz"

	output:
		gff = "database/{id}_genomic.gff"

	shell:
		"gzip -f -d -k {input.gff} > {output.gff}"


#rule to create a database for a specie
rule makedb:
	params:
		id = "{id}"
	input:
		"database/{id}_protein.faa"
	output:
		"database/{id}_protein.faa.pin"
	shell:
		"makeblastdb -in {input} -dbtype prot > {output}"

#run blast against the new database
rule blastp :
	params:
		id = "{id}"
	input:
		index = "database/{id}_protein.faa.pin",
		faa = "database/{id}_protein.faa",
		model = "PGC_DShi.fsa"
	output:
		"blasts/blast{id}.txt"
	shell:
		"blastp -db {input.faa} -query {input.model} -outfmt \"6 std qcovs\" > {output}"


#rule to create a file with the best_hits using blast_output_manager python script
rule manage_blast :
	params:
		id = "{id}"
	input:
		"blasts/blast{id}.txt",
		"table2.csv",
		"database/{id}_genomic.gff"
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
		"database/{id}_protein.faa"
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
        expand("database/{sample}.faa", sample=SAMPLES)
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

	
