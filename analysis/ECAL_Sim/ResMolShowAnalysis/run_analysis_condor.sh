#!/bin/bash

ORIGIN=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL
source $ORIGIN/init_ilcsoft_v02-02-03.sh
echo $PWD                                                                                                                                                                  

FOLDER=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/analysis/ECAL_Sim/ResMolShowAnalysis
if [ ! -e ${FOLDER}/results_folder ]; then
    mkdir ${FOLDER}/results_folder
fi

cd ${FOLDER}/results_folder

#FOLDER=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/analysis/ECAL_Sim/ResMolShowAnalysis

echo "Submit --- > analysis "
condor_submit ${FOLDER}/run_analysis.sub

cd -
sleep 10s



