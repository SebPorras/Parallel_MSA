#!/bin/bash
#SBATCH --job-name=msaS_300
#SBATCH --partition=coursework
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-03:00:00
#SBATCH --output=%x_%A.out  # %x will be replaced with the job name, and %A with the job number


time ./msa ../data/globin/300_seqs_globin 
