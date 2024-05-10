#!/bin/bash

dir=PID_6_60_100_GeV_e_pi_mu

if [ ! -e $dir ]; then
    mkdir $dir
fi

cp analysis_template.C $dir/analysis.C

cd $dir
echo $PWD
sed -i "s/XENERGIESX/3/g" analysis.C
sed -i "s/XENERGYSTRINGX/6., 60., 100/g" analysis.C

for particle in e- pi- mu-
do
    root -q analysis.C\(\"$particle\"\)
done

cd -
