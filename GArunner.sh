#!/bin/bash
chmod -R 755 pyGAsubmit.sh
for graphfilename in *.gpickle; do
	chmod -R 755 $graphfilename;
	for datfilename in *.bin; do
		chmod -R 755 $datfilename;
		sbatch pyGAsubmit.sh $graphfilename $datfilename
	done
	sbatch pyGAsubmit.sh $graphfilename
done