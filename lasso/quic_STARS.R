library(pulsar)
library(huge)
library(data.table)
args <- commandArgs(TRUE)
eachtiss <- args[1]

library(QUIC)
quicr <- function(data, lambda) {
    S    <- cov(data)
    est  <- QUIC::QUIC(S, rho=1, path=lambda, msg=0, tol=1e-2)
    path <-  lapply(seq(length(lambda)), function(i) {
                tmp <- est$X[,,i]; diag(tmp) <- 0
                as(tmp!=0, "lgCMatrix")
    })
    est$path <- path
    est
}

gencode <- fread("/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt", sep = "\t")
genes <- t(gencode[,1]) #id
names <- t(gencode[,6]) #name
names(names) <- genes

modulesDir <- "/work-zfs/abattle4/amy/autism/lasso/pc-correction/expr_corrected_modules/"
networksDir <- "/work-zfs/abattle4/amy/autism/lasso/networks/quic/"
paramsDir <- "/work-zfs/abattle4/amy/autism/lasso/stars-params/"

#normalize reads
normalize <- function(dat) {
        n = nrow(dat)
        p = ncol(dat)
        rank.dat =  dat # matrix for ranking
        for (i in 1:p){
                rank.dat[,i] = rank(dat[,i])
         }
        U = rank.dat/(n+1)
        norm.dat = qnorm(U)
        norm.dat
}

corner_text <- function(text, location="topright"){
	legend(location,legend=text, bty ="n", pch=NA) 
}

print(paste("Creating networks for", eachtiss))

reads <- readRDS(paste(modulesDir, eachtiss, ".RDS", sep = ""))
module.names <- names(reads)
#module.names <- c("20")


pdf(paste(paramsDir, eachtiss, ".Edges.pdf", sep = ""))
module.networks <- lapply(module.names, function(modu, reads, tiss) {
	norm_Reads <- normalize(reads[[modu]]) #SAMPLES X GENES IN MODU
	genes <- names[colnames(norm_Reads)]
        print(paste("At module", modu))

	lmax <- getMaxCov(norm_Reads)
	lams <- getLamPath(lmax, lmax*.05, len=40)
	quicargs  <- list(lambda=lams)
	nc    <- if (.Platform$OS.type == 'unix') 2 else 1
	out.p    <- pulsar(norm_Reads, fun=quicr, fargs=quicargs, rep.num=100, criterion=c('stars', 'gcd'), lb.stars=TRUE, ub.stars=TRUE, seed=10010)

	gcd_index <- get.opt.index(out.p, 'gcd')
	if (gcd_index < 0 || gcd_index > 39) {
		print("SKIPPED")
		NA
	} else {
		opt.index(out.p, 'gcd') <- get.opt.index(out.p, 'gcd')

		#Manually setting stars to a low lambda value to test
		#opt.index(out.p, 'stars') <- 39 
		print(out.p)
		stars_lambda <- lams[opt.index(out.p, 'stars')]
		gcd_lambda <- lams[opt.index(out.p, 'gcd')]
		module_size <- ncol(norm_Reads)

		fit.p <- refit(out.p)
		stars_network <- fit.p$refit$stars
		gcd_network <- fit.p$refit$gcd
	
		perm_stars_network <- stars_network
		perm_gcd_network <- gcd_network 
		stars_node_degree <- rowSums(stars_network == 1) 
		gcd_node_degree <- rowSums(gcd_network == 1) 
		
		row.names(stars_network) <- genes
		row.names(gcd_network) <- genes
		set.seed(10)
		row.names(perm_stars_network) <- sample(genes,length(genes),replace=FALSE)
		row.names(perm_gcd_network) <- sample(genes,length(genes),replace=FALSE) 

		hist(gcd_node_degree, col=rgb(1,0,0,0.5), xlim = c(0,70), ylim = c(0,70), main=paste(tiss, 'module', modu, ': gcd vs. stars edges'), xlab='Node Degrees')
		hist(stars_node_degree, col=rgb(0,0,1,0.5), add=T)
		box()

		corner_text(paste("module size:", module_size, "\nstars lambda:", stars_lambda, "\ngcd lambda:", gcd_lambda, "\ngcd index:", opt.index(out.p, 'gcd')))
		networks <- list(as.matrix(stars_network), as.matrix(gcd_network), as.matrix(perm_stars_network), as.matrix(perm_gcd_network))
		names(networks) <- c("stars", "gcd", "perm_stars", "perm_gcd")
		networks
	}
}, reads, eachtiss)

names(module.networks) <- module.names
saveRDS(module.networks, file = paste(networksDir, eachtiss, ".networks.RDS", sep = ""))
dev.off()
