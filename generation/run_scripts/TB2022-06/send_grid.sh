#!/bin/bash

#particle="e-"
#particle="mu-"
for particle in "neutron" "kaon-" "e-" "mu-" "pi-" "gamma"
do
    for energy in 60 100 #1 2 3 4 5 6 7 8 9 10 60 100 
    do
	#source generic_condor.sh $particle $energy $conf "grid_-40-40_"$particle$energy"GeV.mac"
	# sbatch -N1 -t 1-0 -o slurm_out/slurm-%j.out --partition=htc ./generic_condor.sh $particle $energy $conf 
	./generic_condor.sh $particle $energy 
	#break
    done
done


#    sleep 20h
#    mv data/* /lustre/ific.uv.es/prj/gl/abehep.flc/SiWECAL/SiWECAL_p_AHCAL_Sim/.

