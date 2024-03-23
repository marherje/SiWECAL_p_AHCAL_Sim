#!/bin/bash

set -e

if   [ "$#" -eq 1 ]; then

  ODIR="$1"
  IDIR="${TTHBB_ANALYSIS}"/"${ODIR}"/mvaInput_mvaEvent_mvaEventP_boostedHbb_cp_reco/Nominal/combined

elif [ "$#" -eq 2 ]; then

  ODIR="$1"
  IDIR="$2"

else

  printf "\n%s\n" ">>> ERROR -- invalid list of command-line argument(s):"
  printf "%s\n"   "             [1] path to output directory"
  printf "%s\n\n" "             [2] path to input directory (optional)"
  exit
fi

if [ -d "${ODIR}" ]; then

  printf "\n%s\n\n" ">>> ERROR -- target output directory already exists: ${ODIR}"
  exit
fi

if [ ! -d "${IDIR}" ]; then

  printf "\n%s\n\n" ">>> ERROR -- target input directory not found: ${IDIR}"
  exit
fi

mkdir -p "${ODIR}"

for i_sample in signal background; do

  i_trai="${IDIR}"/"${i_sample}"Training.root
  i_test="${IDIR}"/"${i_sample}"Testing.root

  if [ ! -f "${i_trai}" ]; then printf "\n%s\n\n" ">>> ERROR -- target input file not found: ${i_trai}"; break; fi;
  if [ ! -f "${i_test}" ]; then printf "\n%s\n\n" ">>> ERROR -- target input file not found: ${i_test}"; break; fi;

  hadd "${ODIR}"/"${i_sample}"_TrainingPlusTesting.root "${i_trai}" "${i_test}"

  printf "\n%s\n\n" ">>> created PSO input file: ${ODIR}/${i_sample}_TrainingPlusTesting.root"

  unset -v i_trai i_test

done
unset -v i_sample

unset -v ODIR IDIR
