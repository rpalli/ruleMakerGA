
GADIR=/home/rpalli/ruleMakerGA/
OUTDIR=/home/rpalli/data/5-15-17/
mkdir -p $OUTDIR
cp -r $GADIR $OUTDIR
cd $OUTDIR
chmod -R 755 pyNCsubmit.sh
sbatch pyNCsubmit.sh
