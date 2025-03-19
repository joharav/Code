#!/bin/bash
#SBATCH --account=joharav0
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=10gb
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --output=compstat_triplets_output.log

module load julia  # Load Julia if needed
julia -p 8 compstat_triplets.jl  # Run Julia script with 8 processes

