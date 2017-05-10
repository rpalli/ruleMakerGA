
GADIR=/home/rpalli/ruleMakerGA/
OUTDIR=/home/rpalli/data/5-10-17/
mkdir -p $OUTDIR
cp -r $GADIR $OUTDIR
cd $OUTDIR
sbatch pyNCsubmit.sh
sleep 5

grep all the graphs

for each graph, 
	sbatch pyGAsubmit.sh graphfilename