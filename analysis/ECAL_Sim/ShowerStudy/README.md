# Overview
This folder contains the main scripts for the ECAL study.
# Structure
## run_analysis_condor.sh

- `Def`: simple bash macro to run the analysis in condor. Uses run_analysis_template.sh, run_analysis_template.sub
	      and analysis_template.

## analysis_template_gaus

- `Def`: same that anaylsis_template but only using gaus fitting, runs ~20x faster but has less precision.	     

## res_mol_plots.C

- `Def`: selection of plots covering linearity study, shower profile study, and some shower characterization variables
- `Func`: void res_mol_plots(string particle, bool transformed = true, bool save = false)
- `Par`: particle can be "e-", "pi-", "mu-"
- `Par`: transformed true  means that the x-axis uses 1/sqrt(E) instead of E
- `Par`: save true is for saving the plots in .eps files

## single_histo.C

- `Def`: plot one of the TH1F variables of a specific energy, for e-, pi- and mu- simultaneously
- `Func`: single_histo(TString varname, TString energy, bool save = false)
- `Par`: varname is the name of the TH1F histo selected
- `Par`: energy... is energy
- `Par`: save plot in .eps or not
	
## plots_macro.sh

- `Def`: simple bash macro to farm plots while you can focus on smth else

### test_folder
- `Def`: Self-Explanatory title. Do all the testing there.

# Old code:
## Analysis_singlefile.C

- `Def`: reads the build files and prepare a single ROOT file with all the TH1F, fits and average values.	
	      Use for testing. You'll need to search "XENERGIESX" and replace it properly.
- `Func`: void analysis_template (string particle)
- `Par`: particle can be "e-", "pi-", "mu-"

## run_analysis_singlefile

- `Def`: Same thing but works locally with the output of Analysis_singlefile.C, can be used for testing and/or running the analysis locally.
