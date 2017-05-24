
GADIR=/home/rpalli/ruleMakerGA/
OUTDIR=/home/rpalli/data/in_degree_focus_biased_crossover_40gen/
mkdir -p $OUTDIR
cp -r $GADIR $OUTDIR
cd $OUTDIR
chmod -R 755 ruleMakerGA/pyNCsubmit.sh
chmod -R 755 ruleMakerGA/networkConstructor.py
