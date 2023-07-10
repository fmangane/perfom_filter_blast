
import pandas as pd
import gff3_parser

outblast = pd.read_table(snakemake.input[0],sep = "\t")
outblast.columns=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","coverage"]




#Adding columns to the df
outblast["qlength"]=abs(outblast["qend"]-outblast["qstart"])+1
outblast["slength"]=abs(outblast["send"]-outblast["sstart"])+1

# Print the DataFrame
print(outblast.columns)
print(outblast["coverage"])



#filtering by e-value
for i in range(len(outblast)):
    if outblast.loc[i,"evalue"] > 0.000001 :
        outblast = outblast.drop([i])
outblast.reset_index(inplace = True, drop = True)

for i in range(len(outblast)):
    if outblast.loc[i,'coverage'] < 50:
        outblast = outblast.drop([i])
outblast.reset_index(inplace=True, drop=True)

print(outblast)

if outblast.empty == True:
    print("it's empty")
    print(outblast.columns)
    outblast['qseqid'] = ["nan"]
    outblast['sseqid'] = ["nan"]
    outblast['length'] = ["nan"]
    outblast['mismatch'] = ["nan"]
    outblast['gapopen'] = ["nan"]
    outblast['qstart'] = ["nan"]
    outblast['qend'] = ["nan"]
    outblast['sstart'] = ["nan"]
    outblast['send'] = ["nan"]
    outblast['evalue'] = ["nan"]
    outblast['bitscore'] = ["nan"]
    outblast['coverage'] = ["nan"]
    outblast['qlength'] = ["nan"]
    outblast['slength'] = ["nan"]
    outblast["qGene"] = ["nan"]
    outblast["qFunction"]  = ["nan"]
    outblast["qStrand"]= ["nan"]
    outblast["qStart"]= ["nan"]
    outblast["qEnd"]= ["nan"]
    outblast["sGene"] = ["nan"]
    outblast["sFunction"]  = ["nan"]
    outblast["sStrand"]= ["nan"]
    outblast["sStart"]= ["nan"]
    outblast["sEnd"]= ["nan"]
    outblast["sGenome"] = ["nan"]
    outblast["sAssembly"] = ["nan"]


    outblast["Organism"]= snakemake.params[0]
    print(outblast)
    best_hits_csv = outblast.to_csv(snakemake.output[0], index = True)
    print('\nCSV String:\n', best_hits_csv)

else :

    # selectiong only the best hits for each sseqid 
    droplist=[]
    for i in range(len(outblast)):
        a = outblast.loc[i,"sseqid"]
        b = outblast.loc[i,"qseqid"]
        mac= outblast.loc[i,"bitscore"]
        for j in range(len(outblast)):
            if outblast.loc[j,"sseqid"] == a and outblast.loc[j,"qseqid"] != b and outblast.loc[j,"bitscore"] < mac :
                droplist.append(j)
    outblast = outblast.drop(droplist)
    outblast.reset_index(inplace=True, drop=True)



    print(outblast)



    Table2 = pd.read_csv(snakemake.input[1],header=0)


    #getting qseq gene name from table2

    Genelist = []
    Functionlist = []
    qStrandlist = []
    qStartlist = []
    qEndlist = []
    for i in range(len(outblast)):
        a= outblast.loc[i,"qseqid"]
        for j in range(len(Table2)):
            b = Table2.loc[j,"Accession"]
            if a == b:
                Genelist.append(Table2.loc[j,"Gene"])
                Functionlist.append(Table2.loc[j,"Function"])
                qStrandlist.append(Table2.loc[j,"Strand"])
                qEndlist.append(Table2.loc[j,"qEnd"])
                qStartlist.append(Table2.loc[j,"qStart"])


    #adding them to outblast whether it's empty or not

    outblast["qGene"] = Genelist
    outblast["qFunction"] = Functionlist
    outblast["qStrand"]= qStrandlist
    outblast["qStart"]= qStartlist
    outblast["qEnd"]= qEndlist

    print(outblast)


    #getting sseq gene and function from gff file
    #Parsing the gff file
    gff_file = snakemake.input[2]
    #full_data = pd.read_csv(gff_file,header = 0 , sep = "\t")
    just_tabular = gff3_parser.parse_gff3(gff_file, verbose = True, parse_attributes = False)
    full_data = gff3_parser.parse_gff3(gff_file,verbose = False,  parse_attributes=True)
    print(full_data.columns)
    
    #getting the list of gene in the same order than the accession
    sgenelist = []
    for i in range(len(outblast)):
        acc = outblast.loc[i,'sseqid']
        for j in range(len(full_data)):
            if acc == full_data.loc[j,'Name']:
                sgenelist.append(full_data.loc[j,'gene'])
    #getting the list of functions in the same order than the accession
    sfunctionlist = []
    for i in range(len(outblast)):
        acc = outblast.loc[i, 'sseqid']
        for j in range(len(full_data)):
            if acc == full_data.loc[j,'Name']:
                sfunctionlist.append(full_data.loc[j,'product'])
    #getting the list of stands
    sStrandlist = []
    for i in range(len(outblast)):
        acc = outblast.loc[i, 'sseqid']
        for j in range(len(full_data)):
            if acc == full_data.loc[j,'Name']:
                sStrandlist.append(full_data.loc[j,'Strand'])
    sStartlist = []
    for i in range(len(outblast)):
        acc = outblast.loc[i, 'sseqid']
        for j in range(len(full_data)):
            if acc == full_data.loc[j,'Name']:
                sStartlist.append(full_data.loc[j,'Start'])
    sEndlist = []
    for i in range(len(outblast)):
        acc = outblast.loc[i, 'sseqid']
        for j in range(len(full_data)):
            if acc == full_data.loc[j,'Name']:
                sEndlist.append(full_data.loc[j,'End'])
    sGenomelist = []
    if "genome" in full_data.columns:
	    for i in range(len(outblast)):
		acc = outblast.loc[i,'sseqid']
		for j in range(len(full_data)):
		    if acc == full_data.loc[j,'Name']:
		        sGenomelist.append(full_data.loc[j,'genome'])
    else :
    	for i in range(len(outblast)):
    		sGenomelist.append("nan")
    sSeqIdlist = []
    for i in range(len(outblast)):
        acc = outblast.loc[i,'sseqid']
        for j in range(len(full_data)):
            if acc == full_data.loc[j,'Name']:
                sSeqIdlist.append(full_data.loc[j,'Seqid'])


    #adding them to outblast

    outblast["sGene"] = pd.Series(sgenelist)

    outblast["sProduct"] = pd.Series(sfunctionlist)

    outblast["sStrand"]= pd.Series(sStrandlist)

    outblast["sStart"]= pd.Series(sStartlist)

    outblast["sEnd"]= pd.Series(sEndlist)

    outblast["sGenome"] = pd.Series(sGenomelist)

    outblast["sAssembly"] = pd.Series(sSeqIdlist)

    outblast["Organism"]= snakemake.params[0]
    print(outblast)


    
    #saving the file of the best hits
    best_hits_csv = outblast.to_csv(snakemake.output[0], index = True)
    print('\nCSV String:\n', best_hits_csv)
