#!/bin/sh
#SBATCH --time=40:0:0
#SBATCH --mem=60G
#SBATCH --partition=unlimited

Rscript quic_STARS.R $1 > log/quic_$1.out
