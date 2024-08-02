#!/bin/bash

dir=plots_PID_3_categories

if [ ! -e $dir ]; then
    mkdir $dir
fi

#copydir=~/Pictures/plot_sim/PID

#if [ ! -e $copydir ]; then
#    mkdir $copydir
#fi

for production in PID_0.5_4.5_8.5_GeV_e_nue_neutron #PID_6_GeV_e_pi_mu PID_100_GeV_e_pi_mu PID_6_60_100_GeV_e_pi_mu
do
    for variable in NumHits SumEnergy WSumEnergy Radius90 Barycenter_z Barycenter_r Barycenter_x Barycenter_y NumHits_layer_1 NumHits_layer_2 NumHits_layer_3 NumHits_layer_4 NumHits_layer_5 NumHits_layer_6 NumHits_layer_7 NumHits_layer_9 NumHits_layer_10 NumHits_layer_11 NumHits_layer_12 NumHits_layer_13 NumHits_layer_14 NumHits_layer_n_8 Weighte_layer_8 Weighte_layer_n_8 ShowerNhitMaxLayer ShowerNhitStartLayer ShowerNhitEndLayer ShowerNhitStart10Layer ShowerNhitAverage ShowerNhitMax Radius90_layer_0 Radius90_layer_1 Radius90_layer_2 Radius90_layer_3 Radius90_layer_4 Radius90_layer_5 Radius90_layer_6 Radius90_layer_7 Radius90_layer_8 Radius90_layer_9 Radius90_layer_10 Radius90_layer_11 Radius90_layer_12 Radius90_layer_13 Radius90_layer_14 Sume_layer_0 Sume_layer_1 Sume_layer_2 Sume_layer_3 Sume_layer_4 Sume_layer_5 Sume_layer_6 Sume_layer_7 Sume_layer_8 Sume_layer_9 Sume_layer_10 Sume_layer_11 Sume_layer_12 Sume_layer_13 Sume_layer_14 ShowerSumeMaxLayer ShowerSumeStartLayer ShowerSumeEndLayer ShowerSumeStart10Layer ShowerSumeAverage ShowerSumeMax MIP_Likeness
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

