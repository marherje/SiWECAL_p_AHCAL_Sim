#!/bin/bash
(
source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh

export MARLIN_DLL=$PWD/lib/libDigitization.so:$MARLIN_DLL

Marlin $1
)
