#!/bin/bash                                                                                                                                                                                                                                   

ORIGIN=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL
source $ORIGIN/init_ilcsoft_v02-02-03.sh

FOLDER=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/analysis/ECAL_Sim/ResMolShowAnalysis
cp $FOLDER/analysis_template.C analysis_2.C

dir=$PWD
echo ${dir}
sed -i "s/XenergyX/2/g" analysis_2.C 
echo "running neutron 2 GeV"
root -q analysis_2.C\(\"neutron\"\)

