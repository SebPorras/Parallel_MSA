#!/bin/bash
#SBATCH --job-name=floatfast_250_seqs_globin
#SBATCH --partition=cosc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G # memory (MB)

time ./msa_float_fast ./data/globin/215_seqs_globin 