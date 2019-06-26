#install.packages("QUIC", repos = "http://cran.us.r-project.org")
library(QUIC)

tissues <- c("BRNACC")
args <- commandArgs(TRUE)
t <- args[1] 
modulesdir <- "/work-zfs/abattle4/amy/autism/lasso/modules/"
networksdir <- "/work-zfs/abattle4/amy/autism/lasso/networks/" 

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

networks <- lapply(tissues, function(tissue, eachtiss, modulesDir, networksDir) {
	print(paste("Creating networks for", eachtiss))

	reads <- readRDS(paste(modulesDir, eachtiss, ".matrix.RDS", sep = ""))
	module.names <- names(reads)

	pdf(paste(networksDir, eachtiss, ".Density.pdf", sep = ""))
	module.networks <- lapply(module.names, function(modu, read, tiss) {
		norm_Reads <- normalize(read[[modu]])
		print(paste("At module", modu))
		Psequence <- c(seq(0.05, 0.1, by = 0.005), seq(0.2, 0.8, by = 0.1))	
		networks <- lapply(Psequence, function(p, norm_reads, tis, mod) { 
	        	S <- cov(norm_reads)
        		rho <- matrix(p, nrow = nrow(S), ncol = ncol(S))
        		diag(rho) <- 0

        		results <- QUIC(S, rho, path = NULL, tol = 1e-04, msg = 1, maxIter = 1000, X.init = NULL, W.init = NULL)

        		Theta <- results$X

        		node_degree <- rowSums(Theta != 0)
			total <- length(node_degree)
			#hist(node_degree, main = paste("Density for", tis, "module", mod, "lambda", p), col = "blue", breaks = 15)
			unique_degree <- unique(unlist(node_degree))

			unique_proport <- c()
			for (d in unique_degree) {
				unique_proport <- c(unique_proport, (sum(node_degree >= d)/total))
			}

			LM = lm(log10(unique_proport) ~ log10(unique_degree))
			R_squared = summary(LM)$r.squared

			print(paste("R^2 for lambda", p, ":", R_squared))
			#Theta
			
			avg_nodes <- mean(node_degree)
			avg_nodes
		}, norm_Reads, tiss, modu)

		avg <- c()	
		for (p in networks) {
			avg <- c(avg, p[1])
		}
		x <- Psequence
		plot(x, avg, xlab = "lambda", ylab="average degrees", main = paste("Average Density for", tiss, "module", modu))
                lines(x[order(x)], avg[order(x)], col = "brown")

		names(networks) <- Psequence
		networks
	}, reads, eachtiss)
	dev.off()

	names(module.networks) <- module.names
	#saveRDS(module.networks, file = paste(networksDir, eachtiss, ".networks.RDS", sep = ""))

}, t, modulesdir, networksdir)
print("Completely done")
