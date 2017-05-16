chmod -R 755 pyGAsubmit.sh
for graphfilename in *.gpickle; do
	chmod -R 755 graphfilename;
	sbatch pyGAsubmit.sh graphfilename
done