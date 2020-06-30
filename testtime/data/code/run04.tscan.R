#!/bin/bash -l
#SBATCH --partition=lrgmem
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH -A hji7
ml  R/3.6.1
Rscript 04_test_null.R tscan
