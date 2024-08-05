#!/usr/bin/env python3

# Short sdcript to extract conversion factors from the conversion auxiliary output file
# Assumes 15 layers
# Example usage: ./conversion2steering.py [Aux file from conversion]

import ROOT as r, sys

n_layers = 15
means = []

f = r.TFile.Open(sys.argv[1], "OPEN")

for i_layer in range(n_layers):
    h = f.Get("Energy_GeV__layer_"+str(i_layer))
    means.append(str(h.GetFunction("f").GetParameter(1)))

print("Add to LCIO2Build template the conversion factors as:")
print("<parameter name=\"GeV2MIP\" type=\"string\">" + " ".join(means)+ "</parameter>")