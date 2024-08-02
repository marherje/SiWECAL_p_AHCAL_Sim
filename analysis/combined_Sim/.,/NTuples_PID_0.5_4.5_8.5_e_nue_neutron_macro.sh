#!/bin/bash

dir=PID_0.5_4.5_8.5_GeV_e_nue_neutron

basepath=/lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/ECAL/submit_jobs/LCIO2build_folder/LCIO2build_output/
#\/lustre\/ific\.uv\.es\/prj\/gl\/abehep\.flc\/SiWECAL\/Test\_e\-\_neutron\_nue\/

if [ ! -e $dir ]; then
    mkdir $dir
fi

cp analysis_template.C $dir/analysis.C

cd $dir
echo $PWD
sed -i "s/XENERGIESX/3/g" analysis.C
sed -i "s/XENERGYSTRINGX/0.5, 4.5, 8.5/g" analysis.C
sed -i -e "s:XBASEPATHX:$basepath:g" analysis.C

for particle in e- neutron nu_e
do
    root -q analysis.C\(\"$particle\"\)
done

cd -
