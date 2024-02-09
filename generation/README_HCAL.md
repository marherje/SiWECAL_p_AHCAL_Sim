### How to use the DD4hep model ###

## View the detector and scan material
To use the model:

Initialise the model by doing source init.sh

You can check it by geoDisplay -compact compact/TBModel_SPSMay2018.xml

You can check the material by materialScan compact/TBModel_SPSMay2018.xml 0 0 -100 0 0 100 (it will do a scan from -100 to 100 cm in z direction)

## Run small test
To run a small test of simulation:

$ cd run_sim/smalltest
$ ./run_sim.sh

Don't forget to change the output file path (SIM.outputFile) in the ddsim_steering_00.py

More complex settings is set up in run_sim folder to send jobs to HTCondor (the energy and particle type can be changed)
Here you should change the setting:
- DEST_PATH to the path for the slcio simulation files
