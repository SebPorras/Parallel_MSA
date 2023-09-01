#!/bin/bash
#SBATCH --job-name=200_msa
#SBATCH --partition=cosc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G # memory (MB)

time ./msaOpt ../milestone_1/data/globin/200_seqs_globin 