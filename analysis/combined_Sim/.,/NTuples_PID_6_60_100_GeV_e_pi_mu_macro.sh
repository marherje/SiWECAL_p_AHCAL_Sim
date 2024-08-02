#!/bin/bash

dir=PID_6_60_100_GeV_e_pi_mu

basepath = /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/submit_jobs/LCIO2build_folder/LCIO2build_output/

if [ ! -e $dir ]; then
    mkdir $dir
fi

cp analysis_template.C $dir/analysis.C

cd $dir
echo $PWD
sed -i "s/XENERGIESX/3/g" analysis.C
sed -i "s/XENERGYSTRINGX/6., 60., 100/g" analysis.C
sed -i "s/XBASEPATHX/$basepath/g" analysis.C

for particle in e- pi- mu-
do
    root -q analysis.C\(\"$particle\"\)
done

cd -
