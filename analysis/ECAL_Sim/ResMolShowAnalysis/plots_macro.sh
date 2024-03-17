#!/bin/bash

dir=plots_March

if [ ! -e $dir ]; then
    mkdir $dir
fi

if [ ! -e $dir/general ]; then
    mkdir $dir/general
fi

if [ ! -e $dir/histos ]; then
    mkdir $dir/histos
fi

for particle in e- #pi- mu-
do
    root -q res_mol_plots.C\(\"$particle\",true,true\)
    root -q res_mol_plots.C\(\"$particle\",false,true\)
    mv *eps $dir/general/.
    mv *png $dir/general/.
done

for variable in Barycenter_z
do
    for energy in 2 #20 100
    do
	root -q single_histo.C\(\"$varname\",\"$energy\",true\)
	mv *eps $dir/histos/.
	mv *eps $dir/general/.
    done
done

sleep 1s



