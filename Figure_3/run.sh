#!/bin/bash
#SBATCH --job-name=000
#SBATCH --output=./000.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb
#SBATCH --qos=short
#SBATCH --time=6:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=carolina.monzo@csic.es

source ~/.bashrc
conda deactivate
conda activate seurat

Rscript make_many.R --seed 000
