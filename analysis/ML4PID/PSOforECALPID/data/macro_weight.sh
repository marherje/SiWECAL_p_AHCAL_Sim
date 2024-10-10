#!/bin/bash
rm *~

for prod in PID_6_GeV #PID_100_GeV_e_pi_mu PID_6_60_100_GeV_e_pi_mu
do
    cp AddWeightBranch.C $prod/.
    cd $prod
    echo $prod
    for sample in resolution_e-_result resolution_pi-_result resolution_mu-_result resolution_kaon-_result resolution_gamma_result resolution_neutron_result
    do
	root -q AddWeightBranch.C\(\"$sample\"\)
    done
    cd -
done
