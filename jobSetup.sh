
GADIR=/home/rpalli/ruleMakerGA/
OUTDIR=/home/rpalli/data/5-15-17/
chmod -R 755 /home/rpalli/ruleMakerGA/pyNCsubmit.sh
mkdir -p $OUTDIR
cp -r $GADIR $OUTDIR
cd $OUTDIR
sbatch pyNCsubmit.sh
