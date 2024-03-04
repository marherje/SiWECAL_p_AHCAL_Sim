#!/bin/bash

if [ ! -e results_folder ]; then
    mkdir results_folder
fi

for particle in e- pi- mu-
do
    root analysis.C\(\"$particle\"\)
    mv *${particle}*root results_folder/.
done

sleep 10s
mv *root results_folder/.


