Please read both READMEs from both detectors.

# GEOMETRY
All the geometry decription.
- `ECALHCAL_TB2022-06` has the working one that's being used for all the analysis.
- `ECALHCAL_ILD_IDR` has the one with the full-lenght ILD prototype.

# RUN_SCRIPTS
Everything needed to run the simulation by using the geometries.
TB2022-06 has the test beam geometry prepared.
- `TB2022-06`: `generic_condor.sh` in the template file that runs each particle run.
	    `send_grid.sh` is the launch file.
	    `restart_everything.sh` does what its name says, be careful with it.