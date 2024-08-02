#!/bin/bash

copydir=~/Pictures/plot_sim/PID_Performance

if [ ! -e ~/Pictures/plot_sim ]; then
    mkdir ~/Pictures/plot_sim
fi

if [ ! -e $copydir ]; then
    mkdir $copydir
fi

for production in PID_6_60_100_GeV_e_pi_mu_n
do
    dir=$production

    if [ ! -e $dir ]; then
	mkdir $dir
    fi
    
    cp Performace_plots_template_4_cat.C $dir/Performance_plots.C

    cd $dir
    echo $PWD

    if [ $production == PID_6_60_100_GeV_e_pi_mu_n ]; then
        symA="e-"
        symB="#pi-"
        symC="#mu-"
	symD="n "
        partA=e
        partB=#pi
        partC=#mu
	partD=n
	prod="6 + 60 + 100 GeV"
    fi

    sed -i "s/XproductionX/"$production"/g" Performance_plots.C
    sed -i "s/XprodX/$prod/g" Performance_plots.C
    sed -i "s/XsymAX/$symA/g" Performance_plots.C
    sed -i "s/XsymBX/$symB/g" Performance_plots.C
    sed -i "s/XsymCX/$symC/g" Performance_plots.C
    sed -i "s/XsymDX/$symD/g" Performance_plots.C
    sed -i "s/XpartAX/"$partA"/g" Performance_plots.C
    sed -i "s/XpartBX/"$partB"/g" Performance_plots.C
    sed -i "s/XpartCX/"$partC"/g" Performance_plots.C
    sed -i "s/XpartDX/"$partD"/g" Performance_plots.C
    sed -i "s/XsetAX/$setA/g" Performance_plots.C
    sed -i "s/XsetBX/$setB/g" Performance_plots.C
    sed -i "s/XsetCX/$setC/g" Performance_plots.C
    sed -i "s/XsetDX/$setD/g" Performance_plots.C

    #Cambiar bool a true y descomentar el rm  
    root -q Performance_plots.C\(true\)
    #mv *eps $dir/.
    #mv *png $dir/.
    #rm Performance_plots.C
    cd -
    cp -r $dir $copydir/.
done

sleep 1s
