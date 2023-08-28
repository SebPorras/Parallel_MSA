#!/bin/bash
#SBATCH --job-name=msaOpt_double_250
#SBATCH --partition=cosc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G # memory (MB)

time ./msaOpt_double_fast ../milestone_1/data/globin/215_seqs_globin 