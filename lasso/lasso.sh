#!/bin/sh
#SBATCH --time=45:0:0
#SBATCH --mem=50G
#SBATCH --partition=unlimited

Rscript lasso_BIC.R $1 > log/lasso_BIC_$1.out
