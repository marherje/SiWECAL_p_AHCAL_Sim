#!/bin/bash

if [ ! -e preshower_AHCALLCIO2build_folder ]; then
    mkdir preshower_AHCALLCIO2build_folder
fi
if [ ! -e preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_log ]; then
    mkdir preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_log
fi
if [ ! -e preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_steer ]; then
    mkdir preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_steer
fi
if [ ! -e preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_output ]; then
    mkdir preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_output
fi

template_steer=../steering/TB2022-06_CONF6/preshower_AHCALLCIO2Build.xml
template_sh=preshower_AHCALLCIO2build_gen.sh
template_sub=preshower_AHCALLCIO2build_gen.sub


for particle in e- #neutron kaon- e- pi- mu-
do
  for energy in  6 #60 100 #2 4 6 8 10 20 30 40 50 60 70 80 90 100 125 150 175 200
  do
      filename="preshower_AHCALLCIO2build_"$particle"_"$energy"GeV"
      cp $template_steer preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_steer/$filename".xml"
      sed -i -e 's/XPARTICLEX/'$particle'/g' preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_steer/$filename".xml"
      sed -i -e 's/XENERGYX/'$energy'/g' preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_steer/$filename".xml"
      
      cp $template_sh preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_steer/$filename".sh"
      sed -i -e 's/XcondorfileX/'$filename'/g' preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_steer/$filename".sh"

      cp $template_sub preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_steer/$filename".sub"
      sed -i -e 's/XcondorfileX/'$filename'/g' preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_steer/$filename".sub"

      cd preshower_AHCALLCIO2build_folder/preshower_AHCALLCIO2build_steer
      echo "Submit --- > " $filename
      condor_submit $filename".sub"
      cd -
  done
done

