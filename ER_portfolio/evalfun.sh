#!/bin/bash
#SBATCH --account=joharav0
#SBATCH --job-name=evalfun
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4           # 8 total MPI/Julia processes
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=standard
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joharav@umich.edu

# ------------------------------
# Environment setup
# ------------------------------
module purge
module load julia

export JULIA_NUM_THREADS=4             # threads per process (match ntasks-per-node)
mkdir -p logs

# ------------------------------
# Run Julia job
# ------------------------------
echo "Starting evalfun.jl at $(date) on $(hostname)"
echo "Nodes: $SLURM_JOB_NUM_NODES | Tasks: $SLURM_NTASKS | Threads: $JULIA_NUM_THREADS"

# Each node has 4 tasks → total 8 processes
julia --machine-file $SLURM_NODEFILE -p $SLURM_NTASKS evalfun.jl

echo "Finished evalfun.jl at $(date)"
