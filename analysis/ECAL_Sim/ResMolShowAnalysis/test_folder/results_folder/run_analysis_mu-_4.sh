#!/bin/bash                                                                                                                                                                                                                                   

ORIGIN=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL
source $ORIGIN/init_ilcsoft_v02-02-03.sh

FOLDER=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/analysis/ECAL_Sim/ResMolShowAnalysis/test_folder
cp $FOLDER/analysis_template.C analysis_4.C

dir=$PWD
echo ${dir}
sed -i "s/XenergyX/4/g" analysis_4.C 
echo "running mu- 4 GeV"
root -q analysis_4.C\(\"mu-\"\)

