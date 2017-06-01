
GADIR=/home/rpalli/ruleMakerGA/
OUTDIR=/home/rpalli/data/IFNGunsimplified_rewire/
mkdir -p $OUTDIR
cp -r $GADIR $OUTDIR
cd $OUTDIR
chmod -R 755 ruleMakerGA/pyNCsubmit.sh
chmod -R 755 ruleMakerGA/networkConstructor.py

GADIR=/home/rpalli/ruleMakerGA/
OUTDIR=/home/rpalli/data/IFNGunsimplified_noRewire/
mkdir -p $OUTDIR
cp -r $GADIR $OUTDIR
cd $OUTDIR
chmod -R 755 ruleMakerGA/pyNCsubmit.sh
chmod -R 755 ruleMakerGA/networkConstructor.py
