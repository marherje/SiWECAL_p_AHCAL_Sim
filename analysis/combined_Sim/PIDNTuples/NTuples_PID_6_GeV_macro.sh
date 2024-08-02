#!/bin/bash

dir=PID_6_GeV

basepath=/lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/combined/submit_jobs/combined_LCIO2build_folder/combined_LCIO2build_output/

if [ ! -e $dir ]; then
    mkdir $dir
fi

cp analysis_template.C $dir/analysis.C

cd $dir
echo $PWD
sed -i "s/XENERGIESX/1/g" analysis.C
sed -i "s/XENERGYSTRINGX/6./g" analysis.C
sed -i "s:XBASEPATHX:$basepath:g" analysis.C

for particle in e- pi- mu- kaon- neutron
do
    root -q analysis.C\(\"$particle\"\)
done

cd -
