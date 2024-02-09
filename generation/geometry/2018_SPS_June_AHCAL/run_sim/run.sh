#!/bin/bash
#
#==============================================================
# Running shell script in parallel over multiple cores
#==============================================================

###==============================================================
#set necessary variables

particle=$1
energy=$2
ProcId=$3
echo "ProcId ${ProcId}"

physList="QGSP_BERT_HP"
namesuffix=""
#particle="e-"
#energy=10
energy_sigma=0.0
RUNNUMBER=60972

##DESTINATION PATH
##THIS DIRECTORY SHOULD ALREADY EXIST WHEN THE SCRIPT IS SUBMITTED
DEST_PATH="."

#########=========================================
NEVENTS=NEVENTSPERPROCESS

#########=========================================
i=${ProcId}
let firstEventNumber=${NEVENTS}*${i}
let lastEventNumber=${firstEventNumber}+${NEVENTS}-1
echo "i=${i}, firstEventNumber = ${firstEventNumber}, lastEventNumber = ${lastEventNumber}"

LCIOfile=Run${RUNNUMBER}_${particle}_${energy}GeV_${physList}${namesuffix}_Ev${firstEventNumber}-Ev${lastEventNumber}-Proc${ProcId}.slcio
GEARfile=Run${RUNNUMBER}_${particle}_${energy}GeV_${physList}${namesuffix}_Ev${firstEventNumber}-Ev${lastEventNumber}.gear

echo "${LCIOfile}"

#let seed=${RUNNUMBER}+${lastEventNumber}
seed=${RANDOM}

# sleep to relieve load of Mokka DB
#usleep `expr \( \( ${RANDOM} % 40 \) + 1 \) \* 1333000` ## use this
let sleepTime=(${RANDOM}%40+1)*1333000
usleep sleepTime

#echo the path to the directory where your job runs, in case
#you want to recover files
echo "HOST=${HOST}, TMPDIR=${TMPDIR}, physList=${physList}"
echo "Running simulation for ${particle} at ${energy} GeV, job index ${i}."
echo "Take Mokka executable from ${DD4HEP_PATH}"
echo "Take Mokka env variables from ${DD4HEP_ENV_PATH}"
echo "seed=${seed}"
echo "========================================================================="
echo " "
echo " "

###==============================================================
#change to directory where you run the job
#cd $TMPDIR

#=============================================================#
#                                                             #
#             prepare for the GEANT4 steering file             #
#                                                             #
#==============================================================
rm -rf g4.mac
cat <<EOF > g4.mac
/run/verbose 0
/event/verbose 0
/tracking/verbose 0

#/generator/generator gps
/gps/particle ${particle}
/gps/pos/type Beam
/gps/pos/centre -15 -15 -1000. mm
/gps/pos/sigma_x 27.1 mm
/gps/pos/sigma_y 25.8 mm
/gps/time 0 ns
/gps/direction 0 0 1

##/gps/ene/type Beam
/gps/ene/mono ${energy} GeV
/gps/ene/sigma ${energy_sigma} GeV

/run/beamOn ${NEVENTS}
exit


EOF

#========================================================================#
#                                                                        #
#     run Mokka in batch mode, with the "<<!" and  "!"  Korn-Shell tags  #
#                                                                        #
#========================================================================#
echo "Starting DD4HEP..."

source ../init.sh
source ../bin/thislcgeo.sh

ddsim --runType run --enableG4GPS --macroFile g4.mac --compactFile ../compact/TBModel_SPSJune2018.xml --outputFile ${LCIOfile} --physics.list ${physList} --physics.pdgfile /afs/desy.de/project/ilcsoft/sw/x86_64_gcc48_sl6/v01-17-10/DD4hep/v00-16/examples/DDG4/examples/particle.tbl --enableDetailedShowerMode --random.seed ${seed} --physics.rangecut 0.05 --action.calo Geant4ScintillatorCalorimeterAction

#===============================================================#
#                                                               #
#           copy output files                                   #
#                                                               #
#===============================================================#
if [ -e ${LCIOfile} ]
    then
    cp ${LCIOfile} ${DEST_PATH}
        #if file was copied, delete it from this directory
    if [ $? -eq 0 ]
        then
        rm ${LCIOfile}
    fi
else
    echo "File ${LCIOfile} does not exist...."
fi

#===========================================================#
#                                                           #
#      do some additional clean-up                          #
#                                                           #
#===========================================================#
echo "==================================================="
rm g4.mac
echo "Last look: "
ls -lhtr
echo "Cleaning up... end."
