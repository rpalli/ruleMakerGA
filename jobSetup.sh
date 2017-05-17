
GADIR=/home/rpalli/ruleMakerGA/
OUTDIR=/home/rpalli/data/5-17-17/
mkdir -p $OUTDIR
cp -r $GADIR $OUTDIR
cd $OUTDIR
chmod -R 755 ruleMakerGA/pyNCsubmit.sh
chmod -R 755 ruleMakerGA/networkConstructor.py
