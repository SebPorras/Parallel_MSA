#!/bin/bash
#SBATCH --job-name=400_msaP_mean_dist_upgma
#SBATCH --partition=coursework
#SBATCH --nodes=2                   # Number of nodes
#SBATCH --ntasks=2                  # Number of tasks (usually 1 for OpenMP)
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4           # Number of CPU cores to allocate
#SBATCH --time=0-03:03:00
#SBATCH --output=%x_%A.out  # %x will be replaced with the job name, and %A with the job number

module load compiler-rt/latest
module add mkl/latest
module add mpi/openmpi-x86_64

#I would have expected the module loads to add these, but apparently not
export PATH=/opt/local/stow/cuda-11.1/bin:$PATH
export PATH=/usr/lib64/openmpi/bin:$PATH

time mpiexec -n 2 -map-by node -bind-to none ./msaAvx ./data/globin/400_seqs_globin 
