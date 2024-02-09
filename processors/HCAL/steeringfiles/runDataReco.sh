#!/bin/bash

###
# Skript used for reconstruction of Data, not for simulation output
#
# Skript by  Naoki Tsuji (University of Tokyo), adapted by Erik Buhmann (University of Hamburg)
# August 2018
###

##### Insert runlist and their energy you want to reconstruct
##### if necessary amend input and output path

ENERGY=60
RUN_LIST=(61406
61407
61408
61409
61410
61411)

INPUT_PATH="/pnfs/desy.de/calice/tb-cern/native/cernAhcalJune2018/slcio/Pion/n"$ENERGY"GeV/PP/"
AHCAL_PATH="/pnfs/desy.de/calice/tb-cern/native/cernAhcalJune2018/AhcalRaw/Pion/n"$ENERGY"GeV/PP/"
BIF_PATH="/pnfs/desy.de/calice/tb-cern/native/cernAhcalJune2018/BifRaw/Pion/n"$ENERGY"GeV/PP/"
DWC_PATH="/pnfs/desy.de/calice/tb-cern/native/cernAhcalJune2018/DWC/root/"
OUTPUT_PATH="/nfs/dust/ilc/group/flchcal/AHCAL_Testbeam_SPS_June2018/reco_rootfiles/Pion/n"$ENERGY"GeV/PP/"

#######################

# check if output folder exits, if not, create
mkdir -p $OUTPUT_PATH;


for RUN in ${RUN_LIST[@]}
	do
		echo "Run number : $RUN"
		sleep 1

		myfile=$(find "$INPUT_PATH" -name "*$RUN*.slcio")
		ahcalfile=$(find "$AHCAL_PATH" -name "*$RUN*.raw")
		biffile=$(find "$BIF_PATH" -name "*$RUN*.raw")

		if [ $RUN -gt 60000 ] ; then
			RUN_NET=$(($RUN - 60000))
		else
			RUN_NET=$RUN
		fi
		dwcfile=$(find "$DWC_PATH" -name "*$RUN_NET*.root")

		echo $myfile
		echo $ahcalfile
		echo $biffile
		echo $dwcfile
		sleep 1

		if [ -z "$myfile" ]; then
			echo -e "Run $1 have not found!"
			echo 'Will exit in 5 seconds';
			sleep 5;
			break;
		    else
			 echo -e "Run number FOUND! File: $myfile"
		fi
        
 
                # Create Copy of SplitCollection and insert filepath
		cp SplitCollections.xml "SplitCollections_"$ENERGY"tmp.xml"

		sed -i "s#INPUT_SLCIO#$myfile#g" "SplitCollections_"$ENERGY"tmp.xml"
		sed -i "s#INPUT_AHCAL#$ahcalfile#g" "SplitCollections_"$ENERGY"tmp.xml"
		sed -i "s#INPUT_BIF#$biffile#g" "SplitCollections_"$ENERGY"tmp.xml"
		sed -i "s#INPUT_DWC#$dwcfile#g" "SplitCollections_"$ENERGY"tmp.xml"

		fileName1=${myfile##*/}
		fileName2=${fileName1%.*}

		outSplitFileSlcio=$OUTPUT_PATH$"Split_"$fileName2$".slcio"

		echo OutputFile of Split:
		echo $outSplitFileSlcio

		sleep 2

		sed -i "s#OUTPUT_SPLIT_SLCIO#$outSplitFileSlcio#g" "SplitCollections_"$ENERGY"tmp.xml"


                # Create copy of Reconstruction steering file and insert filepaths
		cp Reconstruction_wPS_TC.xml "Reconstruction_wPS_TC_"$ENERGY"tmp.xml"

		fileName1=${myfile##*/}
		fileName2=${fileName1%.*}
		#fileName=$OUTPUT_PATH$fileName2


		outRecoFileSlcio=$OUTPUT_PATH$"Reco_"$fileName2$".slcio"
		outRecoFileRoot=$OUTPUT_PATH$"Reco_"$fileName2$".root"

		echo OutputFiles of Reco:
		echo $outRecoFileSlcio
		echo $outRecoFileRoot

		sleep 2

		sed -i "s#OUTPUT_SPLIT_SLCIO#$outSplitFileSlcio#g" "Reconstruction_wPS_TC_"$ENERGY"tmp.xml"
		sed -i "s#OUTPUT_RECO_SLCIO#$outRecoFileSlcio#g" "Reconstruction_wPS_TC_"$ENERGY"tmp.xml"
		sed -i "s#OUTPUT_RECO_ROOT#$outRecoFileRoot#g" "Reconstruction_wPS_TC_"$ENERGY"tmp.xml"

echo run splitcollections
		./myMarlin "SplitCollections_"$ENERGY"tmp.xml"
echo run reconstruction
		./myMarlin "Reconstruction_wPS_TC_"$ENERGY"tmp.xml"

                rm $outSplitFileSlcio
                rm "SplitCollections_"$ENERGY"tmp.xml"
                rm "Reconstruction_wPS_TC_"$ENERGY"tmp.xml"
	done
