library(MASS)

rm(list = ls())
expDir <- "/work-zfs/abattle4/amy/autism/matrix-eqtl/data/"
tissues <- c("BRNACC", "BRNAMY", "BRNCDT", "BRNCTXB24", "BRNCTXBA9", "BRNHIP", "BRNHYP", "BRNPUT", "BRNSNA")

quantile_normalization <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
   
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
   
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

normalize <- lapply(tissues, function(eachtiss, expressDir)
{
	print(paste("Normalizing ", eachtiss, sep = ""))
	myfile <- paste(expressDir, eachtiss, "/GE.txt", sep = "")
	dat <- read.table(myfile, row.names = 1, header = TRUE)
	names <- c()
	for (n in colnames(dat)) {
		n <- strsplit(n, "[.]")
		names <- c(names, paste(n[[1]][1], n[[1]][2], sep = "-"))
	}
	colnames(dat) <- names
	qnorm_dat <- quantile_normalization(dat)

	output_file <- paste(expressDir, eachtiss, "/normalized_GE.txt", sep = "")
	write.table(qnorm_dat, file = output_file, quote = FALSE, sep = "\t")	
}, expDir)

