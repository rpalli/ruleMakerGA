#!/bin/sh
#SBATCH --partition=standard 
#SBATCH --time=1:00:00 
#SBATCH --mem-per-cpu=1g
#SBATCH --cpus-per-task=24
export OMP_NUM_THREADS=24
#SBATCH --output=nc_out.txt -N 1 -J ncRun
#SBATCH --nodes=1

module load intelpython/2.7.12
srun python networkConstructor.py