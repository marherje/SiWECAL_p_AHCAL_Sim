#!/bin/bash                                                                                                                                                                                                                                   

ORIGIN=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL
source $ORIGIN/init_ilcsoft_v02-02-03.sh

FOLDER=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/analysis/ECAL_Sim/ResMolShowAnalysis/test_folder
cp $FOLDER/analysis_template.C analysis_60.C

dir=$PWD
echo ${dir}
sed -i "s/XenergyX/60/g" analysis_60.C 
echo "running pi- 60 GeV"
root -q analysis_60.C\(\"pi-\"\)

