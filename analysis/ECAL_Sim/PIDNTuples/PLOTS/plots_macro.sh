#!/bin/bash

dir=plots_PID_6_GeV_e_pi_mu

if [ ! -e $dir ]; then
    mkdir $dir
fi
for production in PID_6_GeV_e_pi_mu PID_6_60_150_GeV_e PID_6_60_150_GeV_e_pi_mu
do
    for variable in NumHits SumEnergy WSumEnergy Radius90 Barycenter_z Barycenter_x Barycenter_y NumHits_layer_1 NumHits_layer_2 NumHits_layer_3 NumHits_layer_4 NumHits_layer_5 NumHits_layer_6 NumHits_layer_7 NumHits_layer_9 NumHits_layer_10 NumHits_layer_11 NumHits_layer_12 NumHits_layer_13 NumHits_layer_14 NumHits_layer_n_8 Weight_layer_8 Weight_layer_n_8 ShowerNhitMaxLayer ShowerNhitStartLayer ShowerNhitEndLayer ShowerNhitStart10Layer ShowerNhitAverage ShowerNhitMax
    do
	cp single_histo_template.C single_histo_$production.C

	sed -i "s/XproductionX/"$production"/g" single_histo_$production.C

	if [ $production == PID_6_GeV_e_pi_mu ]; then
	    resA=e-
            resB=pi-
            resC=mu-
	    partA=e-
	    partB=pi-
	    partC=mu-
	    energyA=6
	    energyB=6
	    energyC=6
	fi

	if [ $production == PID_6_60_150_GeV_e ]; then
	    resA=6_e-
            resB=60_e-
            resC=150_e-
            partA=e-
            partB=e-
            partC=e-
            energyA=6
	    energyB=60
	    energyC=150
        fi

	if [ $production == PID_6_60_150_GeV_e_pi_mu ]; then
            resA=e-
            resB=pi-
            resC=mu-
            partA=e-
            partB=pi-
            partC=mu-
            energyA=6,60,150
            energyB=6,60,150
            energyC=6,60,150
        fi

	sed -i "s/XresAX/"$resA"/g" single_histo_$production.C
        sed -i "s/XresBX/"$resB"/g" single_histo_$production.C
        sed -i "s/XresCX/"$resC"/g" single_histo_$production.C
	sed -i "s/XpartAX/"$partA"/g" single_histo_$production.C
        sed -i "s/XpartBX/"$partB"/g" single_histo_$production.C
        sed -i "s/XpartCX/"$partC"/g" single_histo_$production.C
	sed -i "s/XsetAX/"$partA" ("$energyA" GeV)/g" single_histo_$production.C
        sed -i "s/XsetBX/"$partB" ("$energyB" GeV)/g" single_histo_$production.C
        sed -i "s/XsetCX/"$partC" ("$energyC" GeV)/g" single_histo_$production.C

	root -q single_histo_$production.C\(\"$variable\",true\)
	mv *eps $dir/.
	mv *png $dir/.
    done
done

sleep 1s

