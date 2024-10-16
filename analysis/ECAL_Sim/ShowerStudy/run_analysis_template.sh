#!/bin/bash                                                                                                                                                                                                                                   

ORIGIN=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL
source $ORIGIN/init_ilcsoft_v02-02-03.sh

FOLDER=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/analysis/ECAL_Sim/ShowerStudy
cp $FOLDER/analysis_template.C analysis_XENERGYX.C

dir=$PWD
echo ${dir}
sed -i "s/XenergyX/XENERGYX/g" analysis_XENERGYX.C 
echo "running XPARTICLEX XENERGYX GeV"
root -q analysis_XENERGYX.C\(\"XPARTICLEX\"\)

