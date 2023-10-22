#!/bin/bash
#SBATCH --job-name=10_msa
#SBATCH --partition=coursework
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:01:00
#SBATCH --output=./outputSlurm/

time ./msaAvx data/globin/10_seqs_globin 
