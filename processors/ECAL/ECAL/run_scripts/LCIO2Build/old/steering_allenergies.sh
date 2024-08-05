#!/bin/sh

echo "Run in marlin_run_area !!"

RUNNAME="data_calib3"
PROCESSOR="LCIO2Build"
CONF="1"

#PHYSLIST="FTFP_BERT"
PHYSLIST="QGSP_BERT"
echo "Using $PHYSLIST"

#if [ "$1" = "1" ]; then
#  FITMODE="gauss"
#elif [ "$1" = "2" ]; then
#  FITMODE="landau"
#elif [ "$1" = "3" ]; then
#  FITMODE="langauss"
#else
#  echo "Not a valid mode: $1"
#  exit
#fi

# example filename:
#/data_ilc/flc/jimenez/simulations/TB2022/CONF1/ECAL_QGSP_BERT_conf1_e-_1.0GeV_1.slcio

# TB2022-03
#for energy in e-_3GeV mu-_40GeV mu-_0.4GeV mu-_4GeV;
energies="1.0 1.4 1.8 2.2 2.6 3.0 3.4 3.8 4.2 4.6 5.2 5.6 6.0"
for energy in $energies;
#for part_e in mu-_40GeV;
do
  part_e="e-_${energy}GeV"
  for i in `seq 1 20`;
  do
    STEERING_FILE=../../steering/"$PROCESSOR"/ECAL_"$PHYSLIST"_conf"$CONF"_"$part_e"_"$i".xml
    cp ../../steering/templates/LCIO2Build.xml $STEERING_FILE
    nvim -c "%s/input.slcio/\/data_ilc\/flc\/jimenez\/simulations\/TB2022\/CONF${CONF}\/ECAL_${PHYSLIST}_conf${CONF}_${part_e}_${i}.slcio/|wq" $STEERING_FILE
    nvim -c "%s/output.root/\/data_ilc\/flc\/jimenez\/simulations\/TB2022\/CONF${CONF}\/ECAL_${PHYSLIST}_conf${CONF}_${part_e}_${i}_build.root/|wq" $STEERING_FILE
    nvim -c "%s/name=\"MaxRecordNumber\" value=\"1000\"/name=\"MaxRecordNumber\" value=\"5000\"/|wq" $STEERING_FILE
    echo "Wrote $STEERING_FILE"
    #break
  done
  #break
done
