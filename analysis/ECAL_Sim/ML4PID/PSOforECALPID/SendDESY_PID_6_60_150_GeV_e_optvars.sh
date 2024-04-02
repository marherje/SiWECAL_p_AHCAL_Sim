#!/bin/bash

ORIGIN=$PWD
source $ORIGIN/init_ilcsoft_v02-02-03.sh
#echo $PWD

if [ ! -e $ORIGIN/output ]; then
    mkdir $ORIGIN/output
fi

if [ ! -e $ORIGIN/log ]; then
    mkdir $ORIGIN/log
fi

cd $ORIGIN
. /cvmfs/sft.cern.ch/lcg/releases/gcc/8.2.0-3fa06/x86_64-centos7/setup.sh

export PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_96/Python/2.7.16/x86_64-centos7-gcc8-opt/bin:${PATH}
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_96/Python/2.7.16/x86_64-centos7-gcc8-opt/lib:${LD_LIBRARY_PATH}

export CXX=g++
export CC=gcc
source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/bin/thisroot.sh
# Loading newer ROOT:
#/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.02/x86_64-centos7-gcc48-opt/bin/thisroot.sh

python runPSO.py -c config/PID_6_60_150_GeV_e/PID_6_60_150_GeV_e_optvars.txt -o output/PID_6_60_150_GeV_e_optvars > log/Log_PID_6_60_150_GeV_e_optvars.log

#nohup ./sendtoyific.sh >& run_production.out  &
