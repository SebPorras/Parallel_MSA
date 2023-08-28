#!/bin/bash
#SBATCH --job-name=msa_base_230_seqs_globin
#SBATCH --partition=cosc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G # memory (MB)

./msa_base ./data/globin/230_seqs_globin 