#!/bin/bash                                                                                                                                                                                                                                   

ORIGIN=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL
source $ORIGIN/init_ilcsoft_v02-02-03.sh

FOLDER=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/analysis/ECAL_Sim/ResMolShowAnalysis
cp $FOLDER/analysis_template.C analysis_175.C

dir=$PWD
echo ${dir}
sed -i "s/XenergyX/175/g" analysis_175.C 
echo "running kaon- 175 GeV"
root -q analysis_175.C\(\"kaon-\"\)

