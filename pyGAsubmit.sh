#!/bin/sh
#SBATCH --partition=standard 
#SBATCH -a 1-10
#SBATCH -J GArun
#SBATCH -o GA_out_%A_%a.txt
#SBATCH --mem-per-cpu=100MB
#SBATCH -t 10:00:00
#SBATCH -n 1
#SBATCH -c 24
export OMP_NUM_THREADS=24

module load intelpython/2.7.12
srun python GA.py $1 $SLURM_ARRAY_TASK_ID
echo "ran GA"
for fileName in *.pickle; do
	cp $fileName OUTDIR=/home/rpalli/data/5-16-17/$fileName