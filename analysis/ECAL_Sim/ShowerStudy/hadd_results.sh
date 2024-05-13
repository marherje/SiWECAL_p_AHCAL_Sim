#!/bin/bash

dir=results_folder

if [ ! -e $dir ]; then
    mkdir $dir
fi

hadd -f $dir/resolution_e-_result.root $dir/resolution_e-_2_result.root $dir/resolution_e-_4_result.root $dir/resolution_e-_6_result.root $dir/resolution_e-_8_result.root $dir/resolution_e-_10_result.root $dir/resolution_e-_20_result.root $dir/resolution_e-_30_result.root $dir/resolution_e-_40_result.root $dir/resolution_e-_50_result.root $dir/resolution_e-_60_result.root $dir/resolution_e-_70_result.root $dir/resolution_e-_80_result.root $dir/resolution_e-_90_result.root $dir/resolution_e-_100_result.root $dir/resolution_e-_125_result.root $dir/resolution_e-_150_result.root $dir/resolution_e-_175_result.root $dir/resolution_e-_200_result.root

hadd -f $dir/resolution_pi-_result.root $dir/resolution_pi-_2_result.root $dir/resolution_pi-_4_result.root $dir/resolution_pi-_6_result.root $dir/resolution_pi-_8_result.root $dir/resolution_pi-_10_result.root $dir/resolution_pi-_20_result.root $dir/resolution_pi-_30_result.root $dir/resolution_pi-_40_result.root $dir/resolution_pi-_50_result.root $dir/resolution_pi-_60_result.root $dir/resolution_pi-_70_result.root $dir/resolution_pi-_80_result.root $dir/resolution_pi-_90_result.root $dir/resolution_pi-_100_result.root $dir/resolution_pi-_125_result.root $dir/resolution_pi-_150_result.root $dir/resolution_pi-_175_result.root $dir/resolution_pi-_200_result.root

hadd -f $dir/resolution_mu-_result.root $dir/resolution_mu-_2_result.root $dir/resolution_mu-_4_result.root $dir/resolution_mu-_6_result.root $dir/resolution_mu-_8_result.root $dir/resolution_mu-_10_result.root $dir/resolution_mu-_20_result.root $dir/resolution_mu-_30_result.root $dir/resolution_mu-_40_result.root $dir/resolution_mu-_50_result.root $dir/resolution_mu-_60_result.root $dir/resolution_mu-_70_result.root $dir/resolution_mu-_80_result.root $dir/resolution_mu-_90_result.root $dir/resolution_mu-_100_result.root $dir/resolution_mu-_125_result.root $dir/resolution_mu-_150_result.root $dir/resolution_mu-_175_result.root $dir/resolution_mu-_200_result.root


sleep 1s

