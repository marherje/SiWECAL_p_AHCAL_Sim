#!/bin/bash

dir=plots_PID_3_categories

if [ ! -e $dir ]; then
    mkdir $dir
fi

#copydir=~/Pictures/plot_sim/PID

#if [ ! -e $copydir ]; then
#    mkdir $copydir
#fi

for production in PID_6_GeV_e_pi_mu #PID_1_to_10_GeV_e_pi_mu #PID_6_GeV_e_pi_mu #PID_1_to_10_GeV_e_pi_mu #PID_100_GeV_e_pi_mu PID_6_60_100_GeV_e_pi_mu
do
    for variable in MIP_Likeness #hits_max_distance shower_nhit_max_layer shower_nhit_start_layer shower_nhit_end_layer ecal_interaction nhit weighte mol bar_z radius90_layer_4 
    do
	cp single_histo_template_3_cat.C single_histo_3_temp.C

	sed -i "s/XproductionX/"$production"/g" single_histo_3_temp.C
	
	if [ $production == PID_6_GeV_e_pi_mu ]; then
	    resA=e-
            resB=pi-
            resC=mu-
	    partA=e-
	    partB=pi-
	    partC=mu-
	    setA="e- (6 GeV)"
	    setB="#pi- (6 GeV)"
	    setC="#mu- (6 GeV)"
	fi
	if [ $production == PID_100_GeV_e_pi_mu ]; then
            resA=e-
            resB=pi-
            resC=mu-
            partA=e-
            partB=pi-
            partC=mu-
            setA="e- (100 GeV)"
            setB="#pi- (100 GeV)"
            setC="#mu- (100 GeV)"
        fi
	if [ $production == PID_6_60_100_GeV_e_pi_mu ]; then
            resA=e-
            resB=pi-
            resC=mu-
            partA=e-
            partB=pi-
            partC=mu-
            setA="e- (6+60+100 GeV)"
            setB="#pi- (6+60+100 GeV)"
            setC="#mu- (6+60+100 GeV)"
        fi
	if [ $production == PID_0.5_4.5_8.5_GeV_e_nue_neutron ]; then
            resA=e-
            resB=nu_e
            resC=neutron
            partA=e-
            partB=nu_e
            partC=neutron
            setA="e- (0.5+4.5+8.5 GeV)"
            setB="#nu_{e} (0.5+4.5+8.5 GeV)"
            setC="n (0.5+4.5+8.5 GeV)"
        fi
	
	if [ $production == PID_1_to_10_GeV_e_pi_mu ]; then
            resA=e-
            resB=pi-
            resC=mu-
            partA=e-
            partB=pi-
            partC=mu-
            setA="e- (1-10 GeV)"
            setB="#pi- (1-10 GeV)"
            setC="#mu- (1-10 GeV)"
        fi

	sed -i "s/XresAX/"$resA"/g" single_histo_3_temp.C
        sed -i "s/XresBX/"$resB"/g" single_histo_3_temp.C
        sed -i "s/XresCX/"$resC"/g" single_histo_3_temp.C
	sed -i "s/XpartAX/"$partA"/g" single_histo_3_temp.C
        sed -i "s/XpartBX/"$partB"/g" single_histo_3_temp.C
        sed -i "s/XpartCX/"$partC"/g" single_histo_3_temp.C
	sed -i "s/XsetAX/$setA/g" single_histo_3_temp.C
        sed -i "s/XsetBX/$setB/g" single_histo_3_temp.C
        sed -i "s/XsetCX/$setC/g" single_histo_3_temp.C

	root -q single_histo_3_temp.C\(\"$variable\",true\)
	mv *eps $dir/.
	mv *png $dir/.
	rm single_histo_3_temp.C
    done
done

sleep 1s
#cp -r $dir $copydir/.

