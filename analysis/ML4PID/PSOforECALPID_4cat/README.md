# This is directly copy from my original repository for LCFI+ PID
# Link: https://github.com/marherje/PSOforLCFIPlushttps://github.com/marherje/PSOforLCFIPlushttps://github.com/marherje/PSOforLCFIPlus

# ParticleSwarmOptimization (PSO) for LCFI+

Implementation of the particle swarm algorithm to optimize the hyper-parameters and Input Variables of TMVA Classificators.  
The program will create a swarm of Classifiers in the parameter space.  
Each classifier will try to improve its hyperparameter combination.  
Classifiers move in the parameter space using the information of the whole swarm.  
After N iterations, save the best ones to usi it in the Flavour tagging.
Utilizes batch system to parallelize classifier training.
It can potentially optimise the use of variables.  

Structure:
   config -> Configuration files
   data	  -> Data	   previously prepared with MakeNTuples from lcfi+
   output -> Where the particles information are stored
   conf	  -> To move and save the final .conf files.
   log 	  -> Where the PSO is


IMPORTANT REMARK: This is done for the four categories (A,B,C,D) of LCFI+



Tutorial in Tutorial.md   


This is an adaptation/extension of another PSO originally being a 2-class classifier adapted by Andrej Saibel to work at CMS. 
Link to the original repository by A.S.: https://github.com/Andrej-CMS/ParticleSwarmOptimization.git

Thanks Andrej for the much-needed help for this adaptation :)