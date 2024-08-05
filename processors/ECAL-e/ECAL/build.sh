#!/bin/bash

OLD_DIR=$PWD
CPP_DIR=$(dirname "$0")
cd $CPP_DIR

# unset MARLIN_DLL
# source ./init_ilcsoft.sh

mkdir -p build
cd build
cmake -C $ILCSOFT/ILCSoft.cmake -DCMAKE_BUILD_TYPE=Release ..
make install

cd $OLD_DIR
