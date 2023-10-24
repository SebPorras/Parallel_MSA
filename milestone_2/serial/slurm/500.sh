#!/bin/bash
#SBATCH --job-name=msa_serial
#SBATCH --partition=coursework
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-01:00:00

time ./msa ../data/globin/500_seqs_globin 
