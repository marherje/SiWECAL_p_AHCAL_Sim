#!/bin/bash
rm *~

for prod in Both_PID_6_GeV
do
    cp AddWeightBranch.C $prod/.
    cd $prod
    echo $prod
    for sample in resolution_e-_result.root resolution_pi-_result.root resolution_mu-_result.root resolution_neutron_result.root resolution_kaon-_result.root
    do
	root -q AddWeightBranch.C\(\"$sample\"\)
    done
    cd -
done
