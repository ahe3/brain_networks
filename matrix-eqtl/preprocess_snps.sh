#!/bin/sh
#SBATCH --time=20:0:0
#SBATCH --mem=200G
#SBATCH --partition=unlimited

Rscript chr_preprocess_snps.R $1 > log/chr_snps_$1.out
