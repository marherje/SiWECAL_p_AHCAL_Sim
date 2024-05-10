#!/bin/bash

dir=PID_6_60_100_GeV_e_pi_mu_neutron

if [ ! -e $dir ]; then
    mkdir $dir
fi

cp analysis_template.C $dir/analysis.C

cd $dir
echo $PWD
sed -i "s/XENERGIESX/3/g" analysis.C
sed -i "s/XENERGYSTRINGX/6., 60., 100/g" analysis.C

for particle in e- pi- mu- neutron
do
    root -q analysis.C\(\"$particle\"\)
done

cd -
