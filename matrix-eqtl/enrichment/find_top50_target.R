library(data.table)

ti <- "BRNPUT"
gene <- "SYN2"

mat <- fread("../hits/final.cent.trans.txt")
print("Finished central file")
	
mat <- mat[grepl(gene, unlist(mat[,6])), ]
mat <- mat[order(mat[,5]),]

genes <- c()
rows <- c()
for (row in 1:nrow(mat)) {
	gene <- mat[row, 2]
	if (!(gene %in% genes)) {
		genes <- c(genes, gene)
		rows <- c(rows,row)
		if (length(genes) == 50) { break }
	}
}

filename <- paste(ti, ".top50.", gene, ".cent.txt", sep = "")
write.table(mat[rows,], filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

print("Finished writing tables")
