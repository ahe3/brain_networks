#!/bin/sh
#SBATCH --time=20:0:0
#SBATCH --mem=240G
#SBATCH --partition=unlimited

Rscript filter_fdr.R $1 > log/filter_fdr_$1.out

