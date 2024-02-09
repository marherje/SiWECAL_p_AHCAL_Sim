#!/bin/bash

###
# Skript used for reconstruction of Data, not for simulation output
#
# Skript by  Naoki Tsuji (University of Tokyo), adapted by Erik Buhmann (University of Hamburg)
# August 2018
###

##### Insert runlist and their energy you want to reconstruct
##### if necessary amend input and output path

RUN_LIST=(61406)
OVERLAY_LIST=(61407)

INPUT_PATH="PATH_TO_RECO_DIRECTORY/"
OVERLAY_PATH="PATH_TO_RECO_DIRECTORY/"

OUTPUT_PATH="PATH_TO_OUTPUT_DIRECTORY/"

#######################

# check if output folder exits, if not, create
mkdir -p $OUTPUT_PATH;


for((i=0;i<${#RUN_LIST[@]}; i++))
	do
		RUN=${RUN_LIST[i]}
		OVERLAY=${OVERLAY_LIST[i]}
		echo "Run number : $RUN"
		sleep 1

		myfile=$(find "$INPUT_PATH" -name "Reco*$RUN*.slcio")
		overlayfile=$(find "$OVERLAY_PATH" -name "Reco*$OVERLAY*.slcio")

		echo $myfile
		echo $overlayfile
		sleep 1

		if [ -z "$myfile" ]; then
			echo -e "Run $1 have not found!"
			echo 'Will exit in 5 seconds';
			sleep 5;
			break;
		    else
			 echo -e "Run number FOUND! File: $myfile"
		fi
        
 
		cp Overlay.xml "Overlay_"$ENERGY"tmp.xml"

		outfile="Overlay"$RUN"_"$OVERLAY

		outOverlayFileSlcio=$OUTPUT_PATH$outfile$".slcio"
		outOverlayFileRoot=$OUTPUT_PATH$outfile$".root"

		sleep 2

		sed -i "s#INPUT_OVERLAY_SLCIO#$overlayfile#g" "Overlay_"$ENERGY"tmp.xml"
		sed -i "s#OUTPUT_RECO_SLCIO#$myfile#g" "Overlay_"$ENERGY"tmp.xml"
		sed -i "s#OUTPUT_OVERLAY_SLCIO#$outOverlayFileSlcio#g" "Overlay_"$ENERGY"tmp.xml"
		sed -i "s#OUTPUT_OVERLAY_ROOT#$outOverlayFileRoot#g" "Overlay_"$ENERGY"tmp.xml"

echo run overlay
                Marlin "Overlay_"$ENERGY"tmp.xml"

                rm "Overlay_"$ENERGY"tmp.xml"
	done
