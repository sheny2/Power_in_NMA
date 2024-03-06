#!/bin/bash
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --job-name=YSjob_NMA
#SBATCH --partition=volfovskylab

module load R/4.1.1-rhel8
R CMD BATCH Power_DCC.R