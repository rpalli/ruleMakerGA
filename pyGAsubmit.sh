#!/bin/sh
#SBATCH --partition=standard 
#SBATCH -a 1-3
#SBATCH -J GArun
#SBATCH -o GA_out_%A_%a.txt
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -c 24
export OMP_NUM_THREADS=24

module load intelpython/2.7.12
srun python experiments.py $1 $2 $SLURM_ARRAY_TASK_ID
echo "ran GA"

