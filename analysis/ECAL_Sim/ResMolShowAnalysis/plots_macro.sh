#!/bin/bash

dir=plots_May

if [ ! -e $dir ]; then
    mkdir $dir
fi

if [ ! -e $dir/general ]; then
    mkdir $dir/general
fi

for particle in e- #pi- mu-
do
    root -q res_mol_plots.C\(\"$particle\",true,true\)
    root -q res_mol_plots.C\(\"$particle\",false,true\)
    rm shower*InvSq*
    rm nhit*InvSq*
    rm weight*InvSq*
    mv *eps $dir/general/.
    mv *png $dir/general/.
done
sleep 1s

if [ ! -e ~/Pictures/plot_sim ]; then
    mkdir ~/Pictures/plot_sim
fi

if [ -e ~/Pictures/plot_sim/plots_April ]; then
    rm -r ~/Pictures/plot_sim/plots_April
fi

if [ -e ~/Pictures/plot_sim/plots_May ]; then
    rm -r ~/Pictures/plot_sim/plots_May
fi

cp -r ${dir}/general ~/Pictures/plot_sim/.
#cp -r ${dir}/histos ~/Pictures/plot_sim/.

