#!/bin/bash
#SBATCH --job-name=10_seqs_globin
#SBATCH --partition=coursework
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

time ./msa data/globin/10_seqs_globin 