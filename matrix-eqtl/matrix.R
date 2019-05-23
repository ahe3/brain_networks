library(MatrixEQTL)

#!/usr/bin/env Rscript
args = commandArgs(TRUE)
tissue = args[1]
chro = args[2]
print(paste("Starting matrix-eQTL for", tissue, chro))
## Location of the package with the data files.
#base.dir = find.package('MatrixEQTL');
base.dir = ""
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "data/", tissue, "/SNP.", chro, ".txt", sep="");
#snps_location_file_name = paste(base.dir, "data/", tissue, "/snpsloc.txt", sep="");

# Gene expression file name
expression_file_name = paste(base.dir, "data/", tissue, "/GE.", chro, ".txt", sep="");
#gene_location_file_name = paste(base.dir, "data/geneloc.txt", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste(base.dir, "data/", tissue, "/Covariates.txt", sep="");

# Output file name
#output_file_name_cis = paste("data/", tissue, "/", chro, ".cis.txt", sep = "");
output_file_name = paste("data/", tissue, "/", chro, ".trans.txt", sep = "");

# Only associations significant at this level will be saved
#pvOutputThreshold_cis = 0;
pvOutputThreshold = 1;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
#cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);
print("Finished loading genotype data")
## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);
print("Finished loading expression data")
## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}
print("Finished loading covariates data")

## Run the analysis
#snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
#genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
cvrt = cvrt,
output_file_name = output_file_name,
pvOutputThreshold = pvOutputThreshold,
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE, 
pvalue.hist = TRUE,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
#cat('Detected local eQTLs:', '\n');
#show(me$cis$eqtls)
#cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)
