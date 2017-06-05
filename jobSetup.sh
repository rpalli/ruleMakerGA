
GADIR=/home/rpalli/ruleMakerGA/
OUTDIR=/home/rpalli/data/cut_motifs1/
mkdir -p $OUTDIR
cp -r $GADIR $OUTDIR
cd $OUTDIR
chmod -R 755 ruleMakerGA/pyNCsubmit.sh
chmod -R 755 ruleMakerGA/networkConstructor.py
