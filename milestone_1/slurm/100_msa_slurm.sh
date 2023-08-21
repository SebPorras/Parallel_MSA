#!/bin/bash
#SBATCH --job-name=100_seqs_globin
#SBATCH --partition=coursework
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

time ./msa ./data/globin/100_seqs_globin 