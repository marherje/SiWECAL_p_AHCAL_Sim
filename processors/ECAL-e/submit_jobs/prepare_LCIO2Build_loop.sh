#!/bin/bash

if [ ! -e LCIO2build_folder ]; then
    mkdir LCIO2build_folder
fi
if [ ! -e LCIO2build_folder/LCIO2build_log ]; then
    mkdir LCIO2build_folder/LCIO2build_log
fi
if [ ! -e LCIO2build_folder/LCIO2build_steer ]; then
    mkdir LCIO2build_folder/LCIO2build_steer
fi
if [ ! -e LCIO2build_folder/LCIO2build_output ]; then
    mkdir LCIO2build_folder/LCIO2build_output
fi

template_steer=../steering/LCIO2Build.xml
template_sh=LCIO2build_gen.sh
template_sub=LCIO2build_gen.sub


for particle in e- neutron pi- mu-
do
  for energy in 10 #100 125 150 175 200 #2 4 6 8 10 20 30 40 50 60 70 80 90 100 125 150 175 200
  do
      filename="LCIO2build_"$particle"_"$energy"GeV"
      cp $template_steer LCIO2build_folder/LCIO2build_steer/$filename".xml"
      sed -i -e 's/XPARTICLEX/'$particle'/g' LCIO2build_folder/LCIO2build_steer/$filename".xml"
      sed -i -e 's/XENERGYX/'$energy'/g' LCIO2build_folder/LCIO2build_steer/$filename".xml"
      
      cp $template_sh LCIO2build_folder/LCIO2build_steer/$filename".sh"
      sed -i -e 's/XcondorfileX/'$filename'/g' LCIO2build_folder/LCIO2build_steer/$filename".sh"

      cp $template_sub LCIO2build_folder/LCIO2build_steer/$filename".sub"
      sed -i -e 's/XcondorfileX/'$filename'/g' LCIO2build_folder/LCIO2build_steer/$filename".sub"

      cd LCIO2build_folder/LCIO2build_steer
      echo "Submit --- > " $filename
      condor_submit $filename".sub"
      cd -
  done
done

