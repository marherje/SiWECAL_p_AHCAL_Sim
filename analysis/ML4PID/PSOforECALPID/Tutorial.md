
Tutorial for Flavour tagging in LCFI+ (and prerequisites for PSO):
-Three main algorithms: TrackNTuple,  MakeNTuple, TrainMVA.
       -TrackNTuple: Prepare vertexing info, you probably won't need to re-do it.
       -MakeNTuple: Prepare the data for training. Here's where you use for new MC data for new weights.
       -TrainMVA: Use the data to train the BDTs and prepare the weights you use for FT. This is what we optimize with PSO.
-Pre-requisites:
	-You need "pure" samples for each class (b-jets, c-jets, uds-jets).
	-You can use a processor before MakeNTuple to extract the desired flavour in each case.
	     -Here's such processor (https://github.com/marherje/LCIO_Extraction.git).
	     -Notice the ISR cut values that define what we consider signal: we used k_isr<35GeV for the 250GeV and k_isr<50GeV for the 500GeV samples.
	-After you have the MakeNTuples for b,c and uds you're ready to go to the TrainMVA optimisation using the PSO!

Tutorial to use Particle Swarm Optimization (PSO) for LCFI+:

1) Clone the repository

2) Configure/Adapt: The python interface manages the communication between the particles and is run on a portal (use screen because of long runtimes).
The Training of the BDTs is done on the batch system and is implemented in Particle.C. This file is also recompiled when you start the PSO now.
    2.1) Remove .pyc files in case you recompile the PSO after any change!
    2.2) This optimization requires Python2.7 & C & C++ & ROOT. The executable SendIFIC250_dEdx_G05_allvars_catA.sh 
    	 is an example of how to run at IFIC, it may need to be adapted.

3) Build a conf file to suit your needs. Check examples in /config/
     play around with the swarm parameters (at least 25 particles recommended)
     choose the Figure of Merit to optimize (at the moment ROC Integral Average)
     fix Kolmogorov-Smirnoff and Anderson-Darling cut values
     specify TMVA factory and method options
     declare coordinate space you want to search
     specify trees and files
     Variables the swarm starts with
     //pool of additional Variables the swarm will try
     3.1) Recommendation: Run one iteration with KS & AD cuts at 0 and then check 
     +If the best configuration have any hyperparameter too close to the limits.
     +Which values of KS & AD ensure safe results (ROC values very close in both Train/Test sets)

4) Start the Optimization with
    python RunPSO.py Example_PSOConfig.txt
    5.1) A better way is using executables, like:
    nohup ./SendIFIC250_dEdx_G05_allvars_catA.sh > nohup_250_dEdx_G05_allvars_catA.log
    
    This way is more stable and straight-forward. Adapt it to yourself!

5) After each iteration the ten best classifiers are writen to PSOResult.txt
   The best classifier and all necessary information is written to a .conf file
