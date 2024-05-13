#!/bin/bash

ORIGIN=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL
source $ORIGIN/init_ilcsoft_v02-02-03.sh
echo $PWD                                                                                                                                                                  

FOLDER=$PWD
if [ ! -e ${FOLDER}/results_folder ]; then
    mkdir ${FOLDER}/results_folder
fi

cd ${FOLDER}/results_folder
#FOLDER=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/analysis/ECAL_Sim/ResMolShowAnalysis
for energy in 2 4 6 8 10 20 30 40 50 60 70 80 90 100 125 150 175 200
do
    for particle in e- pi- mu- neutron kaon-
    do
	echo "Submit --- > analysis "$particle $energy " GeV"
	cp ${FOLDER}/run_analysis_template.sh run_analysis_${particle}_${energy}.sh
	cp ${FOLDER}/run_analysis_template.sub run_analysis_${particle}_${energy}.sub

	sed -i "s/XENERGYX/"$energy"/g" run_analysis_${particle}_${energy}.sh
	sed -i "s/XENERGYX/"$energy"/g" run_analysis_${particle}_${energy}.sub
	sed -i "s/XPARTICLEX/"$particle"/g" run_analysis_${particle}_${energy}.sh
	sed -i "s/XPARTICLEX/"$particle"/g" run_analysis_${particle}_${energy}.sub
	condor_submit run_analysis_${particle}_${energy}.sub
    done
done
cd -
sleep 10s



