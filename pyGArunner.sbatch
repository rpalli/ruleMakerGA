#!/bin/sh
#SBATCH --partition=debug 

#SBATCH --time=1:00:00 
#SBATCH --mem-per-cpu=1g
#SBATCH --cpus-per-task=24
export OMP_NUM_THREADS=24
#SBATCH --output=slurm_out.txt -N 1 -J ruleMakerRun
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --array=1-10

GADIR="/home/rpalli/rulemakerGA"
OUTDIR="/public/rpalli/data/5-10-17"

if [! -d $OUTDIR  ]
	mkdir -p $OUTDIR
cp $GADIR $OUTDIR
cd $OUTDIR
module load "intelpython/2.7.12"

srun python 



python GA.py $SLURM_ARRAY_TASK_ID > output.txt
