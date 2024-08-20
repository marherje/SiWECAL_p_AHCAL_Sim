#!/bin/bash

folder=combined_LCIO2build_folder/combined_LCIO2build_output

for particle in e- neutron kaon- pi- mu- gamma
do
  for energy in 60 100 #1 2 3 4 5 6 7 8 9 10 #60 100 
  do
      for detector in LCIO2Build AHCALLCIO2Build
      do
	  root -q check_events.C\(\"$folder/output_${detector}_TB2022-06_${particle}_${energy}GeV.root\",\"$detector\",60000\)
      done
  done
done

