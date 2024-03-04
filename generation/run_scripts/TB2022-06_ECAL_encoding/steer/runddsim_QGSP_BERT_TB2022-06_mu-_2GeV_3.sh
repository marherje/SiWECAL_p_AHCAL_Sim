#!/bin/bash
cp -r /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/generation/run_scripts/TB2022-06_ECAL_encoding/steer/runddsim_QGSP_BERT_TB2022-06_mu-_2GeV_3.* .
source /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/init_ilcsoft_v02-03-02_g103_flags17.sh
ddsim --enableG4GPS --macroFile /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/generation/run_scripts/TB2022-06_ECAL_encoding/macros/grid_-0.0-0.0_mu-_2GeV.mac --steeringFile /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/generation/run_scripts/TB2022-06_ECAL_encoding/steer/runddsim_QGSP_BERT_TB2022-06_mu-_2GeV_3.py
&> /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/generation/run_scripts/TB2022-06_ECAL_encoding/log/QGSP_BERT_TB2022-06_mu-_2GeV_3.log
#tar czvf /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/generation/run_scripts/TB2022-06_ECAL_encoding/TB2022-06_QGSP_BERT_TB2022-06_mu-_2GeV_3.slcio.tar.gz TB2022-06_QGSP_BERT_TB2022-06_mu-_2GeV_3.slcio 
#rm /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/generation/run_scripts/TB2022-06_ECAL_encoding/log/errors_runddsim_QGSP_BERT_TB2022-06_mu-_2GeV_3* /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/generation/run_scripts/TB2022-06_ECAL_encoding/log/outfile_runddsim_QGSP_BERT_TB2022-06_mu-_2GeV_3*
mv steer_path/*runddsim_QGSP_BERT_TB2022-06_mu-_2GeV_3.txt log_path/.
