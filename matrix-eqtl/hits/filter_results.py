import sys

t = sys.argv[1]
#tissues = ["BRNACC", "BRNAMY", "BRNCDT", "BRNCTXB24", "BRNCTXBA9", "BRNHIP", "BRNHYP", "BRNPUT", "BRNSNA"]
tissues = ["BRNACC"]
chros = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

gencode = open("/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt", "r")
next(gencode)
gene_name = {}
for line in gencode:
	line = line.split('\t')
	gene = line[0]
	name = line[5]
	gene_name[gene] = name
gencode.close()
print("finished gencode")

cent_folder = "/work-zfs/abattle4/parsana/autism_networks/data/v8/module_central_genes/"
comm_folder = "/work-zfs/abattle4/parsana/autism_networks/data/v8/module_communicator_genes/"
mod_folder = "../wgcna/modules/"

for ti in tissues:
	print(t)
    
	cent_tss_file = open(cent_folder + t + ".tss.txt", "r")
	comm_tss_file = open(comm_folder + t + ".tss.txt", "r")
	cent_tss = {}
    	comm_tss = {}

    	for line in cent_tss_file:
        	line = line.rstrip()
        	line = line.split('\t')
        	if line[3] not in cent_tss:
        		cent_tss[line[3]] = set()
        	cent_tss[line[3]].add((line[4], line[0], line[2])) #chr: pos, central gene, central gene module
	for line in comm_tss_file:
    		line = line.rstrip()
    		line = line.split('\t')
    		if line[5] not in comm_tss:
            		comm_tss[line[5]] = set()
        	comm_tss[line[5]].add((line[6], line[0], line[4], line[3])) #chr: pos, comm gene, pc module, comm gene module
    	cent_tss_file.close()
    	comm_tss_file.close()

	mod_file = open(mod_folder + t + ".modules.mapped.txt", "r")
    	next(mod_file)
    	gene_mod = {}
    	for line in mod_file:
    		line = line.rstrip().split('\t')
        	gene_mod[line[0]] = line[1]
    	mod_file.close()

    	cent_file = open("all.cent.trans.txt", "w")
    	comm_file = open("all.comm.trans.txt", "w")

    	cent_file.write("SNP\tgene\tbeta\tt-stat\tp-value\tcent_genes\n")
    	comm_file.write("SNP\tgene\tbeta\tt-stat\tp-value\tcomm_genes\n")

    	for chro in chros:
        	print(chro)
        	trans_file = open("../matrix-eqtl/data/" + t + "/" + chro + ".trans.txt", 'r')
        	next(trans_file)
        	for line in trans_file:
        		line = line.rstrip().split("\t")
            		chrom = line[0].split("_")[0]
    			pos = float(line[0].split("_")[1])
    			gene = line[1]

    			cent_genes = set()
            		comm_genes = set()

            		if chrom in cent_tss:
            			for pair in cent_tss[chrom]:
                    			pair_pos = float(pair[0])
                        		pair_gene = pair[1] 
					pair_mod = pair[2] #central gene's module

					target_gene_mod = gene_mod[gene]

                        		if pair_mod == target_gene_mod and abs(pair_pos - pos) <= 100000:
                        			cent_genes.add(gene_name[pair_gene])

            		if chrom in comm_tss:
                		for pair in comm_tss[chrom]:
                			pair_pos = float(pair[0])
                    			pair_gene = pair[1] #communicator gene
					pair_mod = pair[2] #pc module 

					target_gene_mod = gene_mod[gene]

                  			if pair_mod == target_gene_mod and abs(pair_pos - pos) <= 100000:
                    				comm_genes.add(gene_name[pair_gene])

    			if len(cent_genes) != 0:
    				cent_file.write("\t".join(line[0:5]))
    				cent_file.write("\t" + str(cent_genes) + "\n")
    			if len(comm_genes) != 0:
    				comm_file.write("\t".join(line[0:5]))
    				comm_file.write("\t" + str(comm_genes) + "\n")


        	trans_file.close()
    
	cent_file.close()
    	comm_file.close()

print("completely done")
