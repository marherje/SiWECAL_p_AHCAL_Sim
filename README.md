# Overview
This repository is a junction of SiWECAL and AHCAL simulations.
The starting point is the SiWECAL-Sim repository in https://github.com/marherje/SiWECAL-Sim but the structure has changed quite a lot to also leave space for the AHCAL analysis.

# Structure
-'Generation': dd4hep sim of a combined geometry for both detectors.
-'Processors': Marlin processors for MIP conversion, digitization, prepare ROOT files from the .slcio, etc.
-'Masking': For ECAL masking of channels.
-'Analysis': Anaylisis of the build files. Shower studies and PID.

## GENERATION
-'geometry': Geometry .xml files. 
Includes original ECAL and AHCAL separate geometries of old simulations.
Includes the total geometry of 2022-06 TB with SiWECAL+AHCAL.
Includes a ILD-IDR scheme with an ECAL (30 layers) and AHCAL (48 layers) together
Note: The combined geometries are lead by one "compact" file that loads a 50x50x50m world and place the beam and detectors inside.
Beware of redefinitions inside the detectors geometries and overlapping of the geometries.
-'run_scripts': Everything necessary to run the simulations by using the geometries in /geometry/.

## PROCESSORS
-'ECAL processors': Proccesing of from .slcio ECalorimetershits into "real" ECAL data.
        -'HCAL processors': Proccesing of from .slcio HCalorimetershits into "real" HCAL data. (TBD)
	-'Joint processors': (TBD)

## ANALYSIS
-'ECAL_Sim': Shower and PID studies.
            -'ShowerStudy': Construction of variables, shower profile, Molière radius, plots, etc.
            -'PIDNTuples': Optimized code for building the "ttree" for PID studies, includes histograms for all of them and macros to obtain plots.
                         Includes macros for different PID scenarios (3 or 4 particles)
            -'ML4PID': Features a Particle Swarm Optimization (PSO) of hyper-parameters for a BDT-Based particle identification (PID).
                     Requires the NTuples from PIDNTuples.
                     -'ML4PID/PSOforECALPID': 3 categories (3 particles)
                     -'ML4PID/PSOforECALPID_4cat': 4 categories (4 particles)
-'More TBD'

#USEFUL CHECKS:
-'To visualize':
"geoDisplay -compact compactgeometryfile.xml"

-'To check materials, distances and possible overlappings':
"materialScan compactgeometryfile.xml x0 y0 z0 x1 y1 z1"
it will display a list of materials moving in a straight line from (x0,y0,z0) to (x1,y1,z1)