## SVD to compute PC loadings
compute.pc.loadings <- function(dat){
	usv <- svd(scale(dat))
	usv$u
}

pc_correct <- function(dat.expr, est.pc.rm){
	# dat.expr - columns are features/genes, rows are samples
        loadings <- compute.pc.loadings(dat.expr)
        n.pc <- c(1:est.pc.rm)
        print(paste("removing", est.pc.rm, "PCs", nrow(dat.expr)))
        ## use residuals from top n.pc principal components
        dat.expr.adjusted <- lm(dat.expr ~ loadings[,n.pc])$residuals
        dat.expr.adjusted
}

qnormalize <- function(dat){
  n = nrow(dat)
  p = ncol(dat)
  rank.dat =  apply(dat, 2, rank, ties.method = "average") # matrix for ranking
  U = rank.dat/(n+1)
  norm.dat = qnorm(U)
  norm.dat
}

expDir <- "/work-zfs/abattle4/parsana/autism_networks/data/v8/filtered_expression_by_tissue/"
exp.fn <- dir(expDir)
corrDir <- "/work-zfs/abattle4/amy/autism/lasso/pc-correction/expr_corrected/"
tissues <- sapply(exp.fn, function(x) strsplit(x, "\\_expression")[[1]][1])
est.pc.rm <- 10 

moduleDir <- "/work-zfs/abattle4/parsana/autism_networks/data/v8/modules/"
modExDir <- "/work-zfs/abattle4/amy/autism/lasso/pc-correction/expr_corrected_modules/"

module.subsetted <- lapply(tissues, function(eachtiss,expressDir, correctDir, modDir, moduleExprDir)
{
	expr <- read.delim(paste(expressDir, eachtiss, "_expression.txt", sep =""), header = T, row.names = 1)
	corrected_expr <- t(pc_correct(t(expr), 10))
	write.table(corrected_expr, file = paste(correctDir, eachtiss, "_expression_corrected.txt", sep = ""), quote = FALSE, sep = "\t")
		
	module.genes <- read.delim(paste(modDir, eachtiss, ".modules.mapped.txt", sep = ""), stringsAsFactors = F)

        module.genes <- module.genes[-which(module.genes$module == "0"),] # exclude 0 modules

        module.names <- unique(module.genes$module)
        # make a list with each item in the list is the genes present in that module
        module.genes <- lapply(module.names, function(x,y){
                sub.genes <- y$feature[which(y$module == x)]
                sub.genes
                }, module.genes)


        # subset expression data by modules
        dat.expr.sub <- lapply(module.genes, function(x,y){
                y.sub <- y[which(rownames(y) %in% x),]
                y.sub.norm <- qnormalize(t(y.sub))
                rownames(y.sub.norm) <- colnames(y.sub)
                colnames(y.sub.norm) <- rownames(y.sub)
                y.sub.norm
                }, corrected_expr)
        names(dat.expr.sub) <- module.names
        saveRDS(dat.expr.sub, file = paste(moduleExprDir, eachtiss, ".RDS", sep = ""))
        print(paste("Number of modules in", eachtiss, length(dat.expr.sub), "Done"))

}, expDir, corrDir, moduleDir, modExDir)
