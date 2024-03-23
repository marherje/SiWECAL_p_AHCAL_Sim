#!/bin/bash

dir=PID_6_60_150_GeV_e

if [ ! -e $dir ]; then
    mkdir $dir
fi

cp analysis_template.C $dir/analysis.C

cd $dir
echo $PWD
sed -i "s/XENERGIESX/3/g" analysis.C
sed -i "s/XENERGYSTRINGX/6., 60., 150/g" analysis.C

for particle in e- #pi- mu-
do
    root -q analysis.C\(\"$particle\"\)
done

cd -
