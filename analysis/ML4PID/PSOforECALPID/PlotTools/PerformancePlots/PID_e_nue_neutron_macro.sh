#!/bin/bash

#copydir=~/Pictures/plot_sim/PID_Performance

#if [ ! -e ~/Pictures/plot-sim ]; then
#    mkdir ~/Pictures/plot-sim
#fi

#if [ ! -e $copydir ]; then
#    mkdir $copydir
#fi

for production in PID_0.5_4.5_8.5_GeV_e_nue_neutron 
do
    dir=$production

    if [ ! -e $dir ]; then
	mkdir $dir
    fi
    
    cp Performace_plots_template_3_cat.C $dir/Performance_plots.C

    cd $dir
    echo $PWD

    if [ $production == PID_0.5_4.5_8.5_GeV_e_nue_neutron ]; then
	symA="e-"
        symB="#nu_{e}"
        symC="neutron"
        partA=e
        partB=#nu_{e}
        partC=n
        prod="0.5+4.5+8.5 GeV"
    fi

    sed -i "s/XproductionX/"$production"/g" Performance_plots.C
    sed -i "s/XprodX/$prod/g" Performance_plots.C
    sed -i "s/XsymAX/"$symA"/g" Performance_plots.C
    sed -i "s/XsymBX/"$symB"/g" Performance_plots.C
    sed -i "s/XsymCX/"$symC"/g" Performance_plots.C
    sed -i "s/XpartAX/"$partA"/g" Performance_plots.C
    sed -i "s/XpartBX/"$partB"/g" Performance_plots.C
    sed -i "s/XpartCX/"$partC"/g" Performance_plots.C
    sed -i "s/XsetAX/$setA/g" Performance_plots.C
    sed -i "s/XsetBX/$setB/g" Performance_plots.C
    sed -i "s/XsetCX/$setC/g" Performance_plots.C
    
    #Cambiar bool a true y descomentar el rm  
    root -q Performance_plots.C\(true\)
    #mv *eps $dir/.
    #mv *png $dir/.
    #rm Performance_plots.C
    cd -
    #cp -r $dir $copydir/.
done

sleep 1s
