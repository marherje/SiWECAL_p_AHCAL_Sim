#!/bin/bash

dir=PID_100_GeV_e_pi_mu

if [ ! -e $dir ]; then
    mkdir $dir
fi

cp analysis_template.C $dir/analysis.C

cd $dir
echo $PWD
sed -i "s/XENERGIESX/1/g" analysis.C
sed -i "s/XENERGYSTRINGX/100./g" analysis.C

for particle in e- pi- mu-
do
    root -q analysis.C\(\"$particle\"\)
done
rm analysis.C
cd -
