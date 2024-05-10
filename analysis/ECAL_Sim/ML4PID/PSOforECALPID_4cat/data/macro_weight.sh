#!/bin/bash
rm *~

for prod in PID_6_60_100_GeV_e_pi_mu_neutron
do
    cp AddWeightBranch.C $prod/.
    cd $prod
    echo $prod
    for sample in resolution_e-_result.root resolution_pi-_result.root resolution_mu-_result.root resolution_neutron_result.root
    do
	root -q AddWeightBranch.C\(\"$sample\"\)
    done
    cd -
done
