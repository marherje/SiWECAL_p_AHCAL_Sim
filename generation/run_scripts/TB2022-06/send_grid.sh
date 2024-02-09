#!/bin/bash

#particle="e-"
#particle="mu-"
for particle in "e-" "mu-" "pi-"
do
    for energy in 175 200 #10 20 30 40 50 60 70 80 90 100 125 150 175 200
    do
	#source generic_condor.sh $particle $energy $conf "grid_-40-40_"$particle$energy"GeV.mac"
	# sbatch -N1 -t 1-0 -o slurm_out/slurm-%j.out --partition=htc ./generic_condor.sh $particle $energy $conf 
	./generic_condor.sh $particle $energy 
	#break
    done
done
