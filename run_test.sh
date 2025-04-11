#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --job-name=reversible
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=256
#SBATCH -o "test_on_reversible.out"
#SBATCH --nodelist=gpu02

module load matlab

export OMP_NUM_THREADS=256
matlab -batch "test_runner"
