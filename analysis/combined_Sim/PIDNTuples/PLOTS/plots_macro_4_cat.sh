#!/bin/bash

dir=plots_PID_4_categories

if [ ! -e $dir ]; then
    mkdir $dir
fi

#copydir=~/Pictures/plot_sim/PID

#if [ ! -e $copydir ]; then
#    mkdir $copydir
#fi

for production in PID_6_GeV_e_pi_mu_neutron PID_1_to_10_GeV_e_pi_mu #PID_6_60_100_GeV_e_pi_mu_neutron
do
    for variable in total_nhit_layer_20 total_nhit_layer_50 nhit_ratio weighte_ratio total_interaction total_MIP_Likeness ecal_MIP_Likeness hcal_MIP_Likeness total_hits_max_distance total_shower_nhit_max_layer total_shower_nhit_start_layer total_shower_nhit_end_layer ecal_interaction hcal_interaction ecal_nhit ecal_weighte ecal_mol ecal_bar_z ecal_radius90_layer_4 hcal_nhit hcal_weighte hcal_mol hcal_bar_z hcal_radius90_layer_4 total_nhit total_weighte total_mol total_bar_z total_radius90_layer_4
    do
	cp single_histo_template_4_cat.C single_histo_4_temp.C

	sed -i "s/XproductionX/"$production"/g" single_histo_4_temp.C

	if [ $production == PID_6_60_100_GeV_e_pi_mu_neutron ]; then
	    resA=e-
            resB=pi-
            resC=mu-
	    resD=neutron
	    partA=e-
	    partB=pi-
	    partC=mu-
	    partD=neutron
	    setA="e- (6+60+100 GeV)"
	    setB="#pi- (6+60+100 GeV)"
	    setC="#mu- (6+60+100 GeV)"
	    setD="n (6+60+100 GeV)"
	fi

	if [ $production == PID_6_GeV_e_pi_mu_neutron ]; then
            resA=e-
            resB=pi-
            resC=mu-
            resD=neutron
            partA=e-
            partB=pi-
            partC=mu-
            partD=neutron
            setA="e- (6 GeV)"
            setB="#pi- (6 GeV)"
            setC="#mu- (6 GeV)"
            setD="n (6 GeV)"
        fi

	if [ $production == PID_1_to_10_GeV_e_pi_mu_neutron ]; then
            resA=e-
            resB=pi-
            resC=mu-
            resD=neutron
            partA=e-
            partB=pi-
            partC=mu-
            partD=neutron
            setA="e- (1-10 GeV)"
            setB="#pi- (1-10 GeV)"
            setC="#mu- (1-10 GeV)"
            setD="n (1-10 GeV)"
        fi

	sed -i "s/XresAX/"$resA"/g" single_histo_4_temp.C
        sed -i "s/XresBX/"$resB"/g" single_histo_4_temp.C
        sed -i "s/XresCX/"$resC"/g" single_histo_4_temp.C
	sed -i "s/XresDX/"$resD"/g" single_histo_4_temp.C
	sed -i "s/XpartAX/"$partA"/g" single_histo_4_temp.C
        sed -i "s/XpartBX/"$partB"/g" single_histo_4_temp.C
        sed -i "s/XpartCX/"$partC"/g" single_histo_4_temp.C
	sed -i "s/XpartDX/"$partD"/g" single_histo_4_temp.C
	sed -i "s/XsetAX/${setA}/g" single_histo_4_temp.C
        sed -i "s/XsetBX/${setB}/g" single_histo_4_temp.C
        sed -i "s/XsetCX/${setC}/g" single_histo_4_temp.C
	sed -i "s/XsetDX/${setD}/g" single_histo_4_temp.C
	
	root -q single_histo_4_temp.C\(\"$variable\",true\)
	mv *eps $dir/.
	mv *png $dir/.
    done
    rm single_histo_4_temp.C
done

#sleep 1s
#cp -r $dir $copydir/.
