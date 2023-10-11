#!/bin/bash
#SBATCH --job-name=10_msa
#SBATCH --partition=coursework
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:05:00

time ./msaAvx ../milestone_1/data/globin/10_seqs_globin 
