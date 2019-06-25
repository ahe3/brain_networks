#!/bin/sh
#SBATCH --time=90:0:0
#SBATCH --mem=240G
#SBATCH --partition=unlimited

python filter_results.py $1 > log/results_$1.out

