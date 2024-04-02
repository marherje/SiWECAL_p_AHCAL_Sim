#!/bin/bash                                                                                                                                                                                                                                   

ORIGIN=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL
source $ORIGIN/init_ilcsoft_v02-02-03.sh

FOLDER=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/analysis/ECAL_Sim/ResMolShowAnalysis/test_folder
cp $FOLDER/analysis_template.C analysis_90.C

dir=$PWD
echo ${dir}
sed -i "s/XenergyX/90/g" analysis_90.C 
echo "running mu- 90 GeV"
root -q analysis_90.C\(\"mu-\"\)

