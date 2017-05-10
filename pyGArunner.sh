#!/bin/sh
#SBATCH --partition=standard 
#SBATCH --time=1:00:00 
#SBATCH --mem-per-cpu=1g
#SBATCH --cpus-per-task=24
export OMP_NUM_THREADS=24
#SBATCH --output=slurm_out.txt -N 1 -J ruleMakerRun
#SBATCH --nodes=1
echo "working"
mkdir /home/rpalli/data/5-10-17
GADIR=/home/rpalli/rulemakerGA
OUTDIR=/home/rpalli/data/5-10-17

if [! -d $OUTDIR  ]
	mkdir -p $OUTDIR
cp -rf $GADIR $OUTDIR
cd $OUTDIR
module load intelpython/2.7.12
srun python networkConstructor.py
echo "ran python"