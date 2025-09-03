#!/bin/bash
#SBATCH --account=joharav0
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=10gb
#SBATCH --partition=standard
#SBATCH --time=03:00:00
#SBATCH --output=evalfun_output.log

module load julia  # Load Julia if needed

julia -p 8 evalfun.jl  # Run Julia with 8 processes, each with 4 threads

