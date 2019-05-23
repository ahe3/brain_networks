#!/bin/sh
#SBATCH --time=10:0:0
#SBATCH --mem=220G
#SBATCH --partition=unlimited

Rscript matrix.R $1 $2> log/matrix_$1_$2.out
