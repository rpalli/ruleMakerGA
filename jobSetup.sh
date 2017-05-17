
GADIR=/home/rpalli/ruleMakerGA/
OUTDIR=/home/rpalli/data/localsearch100gens/
mkdir -p $OUTDIR
cp -r $GADIR $OUTDIR
cd $OUTDIR
chmod -R 755 ruleMakerGA/pyNCsubmit.sh
chmod -R 755 ruleMakerGA/networkConstructor.py
