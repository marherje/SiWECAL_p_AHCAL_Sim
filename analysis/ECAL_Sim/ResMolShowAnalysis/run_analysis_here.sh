#!/bin/bash

dir=$PWD
echo ${dir}
FOLDER=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/analysis/ECAL_Sim/ResMolShowAnalysis

if [ ! -e ${FOLDER}/results_folder ]; then
    mkdir ${FOLDER}/results_folder
fi

cp ${FOLDER}/analysis.C ${dir}/.

for particle in e- pi- mu-
do
    echo "running analysis for " ${particle} 
    root analysis.C\(\"$particle\"\)
    mv *${particle}*root ${FOLDER}/results_folder/.
done

rm ${dir}/analysis.C
rm ${dir}/marquezh.cc
