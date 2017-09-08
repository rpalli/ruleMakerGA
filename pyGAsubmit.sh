#!/bin/sh
#SBATCH --partition=standard 
#SBATCH -J GArun
#SBATCH -o GA_out_%A_%a.txt
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -c 2
export OMP_NUM_THREADS=24

module load intelpython/2.7.12
srun python GA.py $1
echo "ran GA"

