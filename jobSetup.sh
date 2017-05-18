
GADIR=/home/rpalli/ruleMakerGA/
OUTDIR=/home/rpalli/data/localsearch40genCorrected/
mkdir -p $OUTDIR
cp -r $GADIR $OUTDIR
cd $OUTDIR
chmod -R 755 ruleMakerGA/pyNCsubmit.sh
chmod -R 755 ruleMakerGA/networkConstructor.py
