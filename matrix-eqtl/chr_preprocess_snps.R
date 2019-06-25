library(data.table)
args <- commandArgs(TRUE)
t <- args[1] 

chro <-c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")

folder <- "data/"
samples_folder <- "/work-zfs/abattle4/parsana/autism_networks/data/v8/samples/"
geno_folder <- "/work-zfs/abattle4/parsana/autism_networks/data/v8/genotypes/"

print(paste("Preprocessing snps for" , t))

samples <- c()
samples_file <- file(paste(samples_folder, t, ".samples.txt", sep = ""), open = "r")
lines <- readLines(samples_file)
for (line in lines) {
	line <- strsplit(line, "[.]")
	sample <- paste(line[[1]][1], line[[1]][2], sep = "-")
	samples <- c(samples, sample)
}


samples_geno_file <- file(paste(geno_folder, "chr1.012.indv", sep = ""), open = "r")
all_samples <- readLines(samples_geno_file)


snps_file <- file(paste(folder, t, "/all_snps.txt", sep = ""), open = "r")
snps <- readLines(snps_file)

for (c in chro) { #For each chromosome
	print(paste("Appending", c))
	snps_geno_file <- file(paste(geno_folder, c, ".012.pos", sep = ""), open = "r")
	all_snps <- readLines(snps_geno_file)
	snps_geno <- intersect(all_snps, snps)
	snps <- setdiff(snps, snps_geno)
	
	print("Reading matrix")
	mat <- fread(paste(geno_folder, c, ".012", sep = ""))
	mat <- mat[, !1] 
	colnames(mat) <- all_snps

	print("Subsetting matrix")
	mat <- mat[, ..snps_geno]
	mat <- data.matrix(mat)
	rownames(mat) <- all_samples
	mat <- mat[samples,]
	mat <- t(mat)
	output_file <-  paste(folder, t, "/SNP.", c, ".txt", sep = "")
	write.table(mat, file = output_file, quote = FALSE, sep = "\t")
	mat <- NULL
	print(paste("Finished matrix for", c))
} #End For 
	
print(paste("Finished processing snps for" , t))
