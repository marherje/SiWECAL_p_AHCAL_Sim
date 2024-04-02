#!/bin/bash

dir=plots_April

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
    rm shower*InvSq*
    rm nhit*InvSq*
    rm weight*InvSq*
    mv *eps $dir/general/.
    mv *png $dir/general/.
done

for variable in Barycenter_z Barycenter_x Barycenter_y NumHits_layer_1 NumHits_layer_2 NumHits_layer_3 NumHits_layer_4 NumHits_layer_5 NumHits_layer_6 NumHits_layer_7 NumHits_layer_9 NumHits_layer_10 NumHits_layer_11 NumHits_layer_12 NumHits_layer_13 NumHits_layer_14 NumHits_layer_n_8 Weight_layer_8 Weight_layer_n_8 ShowerNhitMaxLayer ShowerNhitStartLayer ShowerNhitEndLayer ShowerNhitStart10Layer ShowerNhitAverage ShowerNhitMax
do
    for energy in 6 60 150
    do
	root -q single_histo.C\(\"$variable\",\"$energy\",true\)
	mv *eps $dir/histos/.
	mv *png $dir/histos/.
    done
done

sleep 1s

