chro = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

tissues = ["BRNACC", "BRNAMY", "BRNCDT", "BRNCTXB24", "BRNCTXBA9", "BRNHIP", "BRNHYP", "BRNPUT", "BRNSNA"]

folder = "/work-zfs/abattle4/amy/autism/matrix-eqtl/data/"

gencode = open("/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt", "r")
next(gencode)
gene_loc = {}
for line in gencode:
        line = line.split('\t')
        gene = line[0]
        chrom = line[1]
	gene_loc[gene] = chrom
gencode.close()
print("finished gencode")

for t in tissues:
	print("Processing " + t)
	file = open(folder + t + "/normalized_GE.txt", "r")
	myfiles = {}
	for c in chro:
		myfiles[c] = open(folder + t + "/GE." + c + ".txt", "w")
	linenum = 1
	for line in file:
		if linenum == 1:
			linenum = 2
			for c in myfiles:
				myfiles[c].write("id\t" + line)
		else:
			gene = line.split('\t')[0]
			loc = gene_loc[gene]
			for c in myfiles:
				if c != loc:
					myfiles[c].write(line)
	for c in myfiles:
		myfiles[c].close()	
	file.close()
		
	
