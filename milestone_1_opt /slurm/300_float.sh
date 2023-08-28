#!/bin/bash
#SBATCH --job-name=msaOpt_float_fast_300
#SBATCH --partition=cosc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G # memory (MB)

./msaOpt_float_fast ../milestone_1/data/globin/230_seqs_globin 