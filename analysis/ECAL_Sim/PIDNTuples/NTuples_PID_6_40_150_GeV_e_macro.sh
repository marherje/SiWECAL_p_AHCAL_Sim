#!/bin/bash

dir=PID_6_40_150_GeV_e

if [ ! -e $dir ]; then
    mkdir $dir
fi

cd $dir
echo $PWD

for energy in 6 40 150
do
    for particle in e- #pi- mu-
    do
	cp ../analysis_template.C analysis_$energy.C
	sed -i "s/XENERGIESX/1/g" analysis_$energy.C
	sed -i "s/XENERGYSTRINGX/"${energy}"./g" analysis_${energy}.C
	sed -i "s/analysis/analysis_"${energy}"/g" analysis_${energy}.C
	root -q analysis_$energy.C\(\"${particle}\"\)	
	mv resolution_e-_result.root resolution_${energy}_e-_result.root 
	rm analysis_$energy.C
    done
done

cd -
