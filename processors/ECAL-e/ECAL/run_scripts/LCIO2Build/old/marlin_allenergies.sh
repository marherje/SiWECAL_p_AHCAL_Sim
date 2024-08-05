#!/bin/sh

PROCESSOR="LCIO2Build"
CONF="1"

#PHYSLIST="FTFP_BERT"
PHYSLIST="QGSP_BERT"
echo "Using $PHYSLIST"

energies="1.0 1.4 1.8 2.2 2.6 3.0 3.4 3.8 4.2 4.6 5.2 5.6 6.0"
for energy in $energies;
do
  part_e="e-_${energy}GeV"
  for i in `seq 20`;
  do
    ~/Projects/Marlin/bin/Marlin ../../steering/"$PROCESSOR"/ECAL_"$PHYSLIST"_conf"$CONF"_"$part_e"_"$i".xml
    #break
  done
  #break
done
