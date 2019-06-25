#!/bin/sh
#SBATCH --time=15:0:0
#SBATCH --mem=50G
#SBATCH --partition=unlimited

Rscript find_top50_target.R > find_top50_target.out

