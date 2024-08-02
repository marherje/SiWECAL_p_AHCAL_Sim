#!/bin/bash
rm *~

for prod in PID_6_GeV_e_pi_mu PID_100_GeV_e_pi_mu PID_6_60_100_GeV_e_pi_mu
do
    cp AddWeightBranch.C $prod/.
    cd $prod
    echo $prod
    for sample in resolution_e-_result.root resolution_pi-_result.root resolution_mu-_result.root
    do
	root -q AddWeightBranch.C\(\"$sample\"\)
    done
    cd -
done
