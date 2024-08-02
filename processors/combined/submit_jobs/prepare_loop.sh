#!/bin/bash

if [ ! -e combined_LCIO2build_folder ]; then
    mkdir combined_LCIO2build_folder
fi
if [ ! -e combined_LCIO2build_folder/combined_LCIO2build_log ]; then
    mkdir combined_LCIO2build_folder/combined_LCIO2build_log
fi
if [ ! -e combined_LCIO2build_folder/combined_LCIO2build_steer ]; then
    mkdir combined_LCIO2build_folder/combined_LCIO2build_steer
fi
if [ ! -e combined_LCIO2build_folder/combined_LCIO2build_output ]; then
    mkdir combined_LCIO2build_folder/combined_LCIO2build_output
fi

template_steer=../steering/TB2022-06_CONF6/combined_LCIO2Build.xml
template_sh=combined_LCIO2build_gen.sh
template_sub=combined_LCIO2build_gen.sub


for particle in e- neutron kaon- pi- mu-
do
  for energy in 6 60 100 #100 125 150 175 200 #2 4 6 8 10 20 30 40 50 60 70 80 90 100 125 150 175 200
  do
      filename="combined_LCIO2build_"$particle"_"$energy"GeV"
      cp $template_steer combined_LCIO2build_folder/combined_LCIO2build_steer/$filename".xml"
      sed -i -e 's/XPARTICLEX/'$particle'/g' combined_LCIO2build_folder/combined_LCIO2build_steer/$filename".xml"
      sed -i -e 's/XENERGYX/'$energy'/g' combined_LCIO2build_folder/combined_LCIO2build_steer/$filename".xml"
      
      cp $template_sh combined_LCIO2build_folder/combined_LCIO2build_steer/$filename".sh"
      sed -i -e 's/XcondorfileX/'$filename'/g' combined_LCIO2build_folder/combined_LCIO2build_steer/$filename".sh"

      cp $template_sub combined_LCIO2build_folder/combined_LCIO2build_steer/$filename".sub"
      sed -i -e 's/XcondorfileX/'$filename'/g' combined_LCIO2build_folder/combined_LCIO2build_steer/$filename".sub"

      cd combined_LCIO2build_folder/combined_LCIO2build_steer
      echo "Submit --- > " $filename
      condor_submit $filename".sub"
      cd -
  done
done

