#!/bin/bash

ORIGIN=/lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/ECAL-e
source $ORIGIN/init_ilcsoft_v02-02-03.sh
echo $PWD                                                                                                                                                                  

#cd $ORIGIN
#. /cvmfs/sft.cern.ch/lcg/releases/gcc/8.2.0-3fa06/x86_64-centos7/setup.sh

#export PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_96/Python/2.7.16/x86_64-centos7-gcc8-opt/bin:${PATH}
#export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_96/Python/2.7.16/x86_64-centos7-gcc8-opt/lib:${LD_LIBRARY_PATH}

#export CXX=g++
#export CC=gcc
#source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/bin/thisroot.sh

#export MARLIN_DLL=$PWD/lib/libDigitization.so:$MARLIN_DLL
#export MARLIN_DLL=$PWD/lib/libConversionProcessor.so:$MARLIN_DLL
export MARLIN_DLL=$ORIGIN/lib/libLCIO2BuildProcessor.so:$MARLIN_DLL


Marlin $ORIGIN/submit_jobs/LCIO2build_folder/LCIO2build_steer/LCIO2build_neutron_10GeV.xml

mv *root $ORIGIN/submit_jobs/LCIO2build_folder/LCIO2build_output/.
mv outfile* $ORIGIN/submit_jobs/LCIO2build_folder/LCIO2build_log/.
mv *log $ORIGIN/submit_jobs/LCIO2build_folder/LCIO2build_log/.
mv error* $ORIGIN/submit_jobs/LCIO2build_folder/LCIO2build_log/.
