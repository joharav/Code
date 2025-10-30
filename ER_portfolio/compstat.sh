#!/bin/bash
#SBATCH --account=joharav0
#SBATCH --job-name=compstat
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4           # 8 total processes
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joharav@umich.edu

# ------------------------------
# Environment setup
# ------------------------------
module purge
module load julia

export JULIA_NUM_THREADS=4
mkdir -p logs

# ------------------------------
# Run Julia job
# ------------------------------
echo "Starting compstat.jl at $(date) on $(hostname)"
echo "Nodes: $SLURM_JOB_NUM_NODES | Tasks: $SLURM_NTASKS | Threads: $JULIA_NUM_THREADS"

# Use all 8 Julia workers across nodes
julia --machine-file $SLURM_NODEFILE -p $SLURM_NTASKS compstat.jl

echo "Finished compstat.jl at $(date)"
