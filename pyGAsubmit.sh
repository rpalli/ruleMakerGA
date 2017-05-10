#!/bin/sh
#SBATCH --partition=standard 
#SBATCH --time=1:00:00 
#SBATCH --mem-per-cpu=1g
#SBATCH --cpus-per-task=24
export OMP_NUM_THREADS=24
#SBATCH --output=GA_out-%j.txt -N 1 -J ruleMakerRun
#SBATCH --nodes=1
#SBATCH --array=1-10

module load intelpython/2.7.12
srun python GA.py $1 $SLURM_ARRAY_TASK_ID
echo "ran GA"