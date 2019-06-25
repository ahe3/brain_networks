library(data.table)
args <- commandArgs(TRUE)
t <- args[1]

#tissues <- c("BRNACC", "BRNAMY", "BRNCDT", "BRNCTXB24", "BRNCTXBA9", "BRNHIP", "BRNHYP", "BRNPUT", "BRNSNA")
tissues <- c("BRNACC")

new <- TRUE
for (ti in tissues) {
	print(t)

	mat <- fread("all.cent.trans.txt", sep = "")
	print("Finished central file")
	pvals <- as.double(unlist(mat[,5]))
	fdr <- p.adjust(pvals, "BH")
	tissue <- rep_len(t, nrow(mat))
	print("Adjusted central gene fdr")

	new_mat <- cbind(mat[,1:5], fdr, mat[,6], tissue)
	mat <- NULL
	new_mat2 <- new_mat[new_mat$fdr <= 0.1,]
	print(paste("Adding", nrow(new_mat2), "to central file"))
	new_mat <- NULL

	if (new) {
		cent_mat <- new_mat2
	} else {
		cent_mat <- rbind(cent_mat, new_mat2)
	}
	new_mat2 <- NULL

	write.table(cent_mat, file = "final_cent_trans_fdr0.1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
	
	mat <- fread("all.comm.trans.txt", sep = "")
	print("Finished communicator file")
        pvals <- as.double(unlist(mat[,5]))
        fdr <- p.adjust(pvals, "BH")
	tissue <- rep_len(t, nrow(mat))
	print("Adjusted communicator gene fdr")

        new_mat <- cbind(mat[,1:5], fdr, mat[,6], tissue)
	mat <- NULL
        new_mat2 <- new_mat[new_mat$fdr <= 0.1,]
	print(paste("Adding", nrow(new_mat2), "to communicator file"))
	new_mat <- NULL

        if (new) {
                comm_mat <- new_mat2
        } else {
                comm_mat <- rbind(comm_mat, new_mat2)
        }
	new_mat2 <- NULL

	write.table(comm_mat, file = "final_comm_trans_fdr0.1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
	new <- FALSE
}

print("Finished writing tables")
