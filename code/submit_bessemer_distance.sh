#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=5G
#SBATCH --time=40:00:00
#SBATCH --mail-user=s.c.mills@sheffield.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --output=run_info/output_distance_%j.txt

# note: stan has been compiled with this GCC version; important 
# version gets specified to this (rather than defaulting to older
# compiler
module load GCC/9.3.0
module load R/4.0.0-foss-2020a

# chains and threads
nchains=2
nthreads=20

# print run information
echo "nchains: $nchains"
echo "nthreads: $nthreads"

# job
Rscript code/run_model_distance.R $nchains $nthreads
