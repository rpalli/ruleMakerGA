
GADIR=/home/rpalli/ruleMakerGA/
OUTDIR=/home/rpalli/data/5-15-17/
git clone https://github.com/rpalli/ruleMakerGA.git
mkdir -p $OUTDIR
cp -r $GADIR $OUTDIR
cd $OUTDIR
sbatch pyNCsubmit.sh
