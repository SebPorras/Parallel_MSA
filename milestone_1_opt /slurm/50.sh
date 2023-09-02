#!/bin/bash
#SBATCH --job-name=50_msa
#SBATCH --partition=cosc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G # memory (MB)

time ./msa_O3 ../milestone_1/data/globin/50_seqs_globin 