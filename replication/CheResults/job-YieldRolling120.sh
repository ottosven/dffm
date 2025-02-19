#!/bin/bash -l 
#SBATCH --job-name Yield120
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=3gb 
#SBATCH --time=24:00:00 
#SBATCH --output=./CheResults/%x-%j.out 
#SBATCH --account=sotto2 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sotto2@uni-koeln.de

module use -a /opt/rrzk/modules/experimental
module load R/4.2.1_intel_mkl
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.2

# R with my R program with MPI workers on multiple nodes
mpirun -quiet -np 1 R --vanilla -f CheYieldRolling120.R
