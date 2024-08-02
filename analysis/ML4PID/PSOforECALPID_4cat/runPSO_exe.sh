#!/bin/bash

set -e

if [ "$#" -gt 1 ]; then

  ODIR="$1"

else

  printf "\n%s\n" " >>> ERROR -- invalid list of command-line argument(s):"
  printf "%s\n"   "           [1]  path to output directory"
  printf "%s\n\n" "           [2+] PSO configuration file(s)"
  exit
fi

if [ -d "${ODIR}" ]; then

  printf "\n%s\n\n" " >>> ERROR -- target output directory already exists: ${ODIR}"
  exit
fi

###

for cfgfile in "${@:2}"; do

  if [ ! -f "${cfgfile}" ]; then

    printf "\n%s\n\n" " >>> WARNING -- target configuration file not found: ${cfgfile}"
    continue
  fi

  if [ ! -d "${ODIR}" ]; then mkdir -p "${ODIR}"; fi;

  if [[ ${cfgfile} != *_config.txt ]]; then continue; fi;

  cfgname=$(basename "${cfgfile}")
  cfgname="${cfgname%_config.txt}"

  printf "%s\n" "> ${ODIR}/${cfgname}"

  nohup ./runPSO.py \
    -c           "${cfgfile}" \
    -o "${ODIR}"/"${cfgname}" \
    >& "${ODIR}"/"${cfgname}".log &

  unset -v cfgname

done
unset -v cfgfile

###

unset -v ODIR
