#!/bin/bash                                                                                                                                                                                                                                   

ORIGIN=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL
source $ORIGIN/init_ilcsoft_v02-02-03.sh

FOLDER=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/analysis/ECAL_Sim/ResMolShowAnalysis
cp $FOLDER/analysis_template.C analysis_70.C

dir=$PWD
echo ${dir}
sed -i "s/XenergyX/70/g" analysis_70.C 
echo "running kaon- 70 GeV"
root -q analysis_70.C\(\"kaon-\"\)

