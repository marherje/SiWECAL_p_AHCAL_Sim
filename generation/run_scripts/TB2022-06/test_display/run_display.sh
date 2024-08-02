#!/bin/bash

ORIGIN=$PWD
source $ORIGIN/init_ilcsoft_v02-03-03.sh

export LD_LIBRARY_PATH=$ORIGIN/lib:$LD_LIBRARY_PATH

glced &

Marlin display.xml
