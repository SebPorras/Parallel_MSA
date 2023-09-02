#!/bin/bash
#SBATCH --job-name=230_msa
#SBATCH --partition=cosc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G # memory (MB)

time ./msa_roll ../milestone_1/data/globin/230_seqs_globin 