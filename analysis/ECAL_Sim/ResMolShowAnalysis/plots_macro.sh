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

for particle in e- pi- mu-
do
    root -q res_mol_plots.C\(\"$particle\",true,true\)
    root -q res_mol_plots.C\(\"$particle\",false,true\)
    mv *eps $dir/general/.
    mv *png $dir/general/.
done

for variable in Barycenter_z Barycenter_x Barycenter_y NumHits_layer_8 NumHits_layer_n_8 Weight_layer_8 Weight_layer_n_8 ShowerNhitMaxLayer ShowerNhitStartLayer ShowerNhitEndLayer ShowerNhitStart10Layer ShowerNhitAverage ShowerNhitMax
do
    for energy in 2 6 20 40 100 150
    do
	root -q single_histo.C\(\"$variable\",\"$energy\",true\)
	mv *eps $dir/histos/.
	mv *png $dir/histos/.
    done
done

sleep 1s

