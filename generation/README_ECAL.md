# Generation of beam events and the SiW-ECAL

The generation is based on the [DD4hep](https://github.com/iLCSoft/lcgeo) models for detectors.

## Setup

The software here has been tested with LCIO: `source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh`

## Structure

- `geometry`: Contains the definition of the detector model geometries and configurations.
- `run_scripts`: A set of tools for running the simulations with different beam conditions. Includes scripts for sending batch jobs.

## Test beam configurations

A script to write the geometry description for a given TB year and configuration
```bash
./write_geom_desc.py TB2022-03_CONF1 > geometry/TB2022-03/ECAL_CONF1.xml 
```
### Configurations

Defined in `confs.py`. A configuration has per layer tree entries: [Tungsten thickness[mm], Seensitive layer (bool), Silicon thickness [µm]]. Per test beam session:

- 2017-06 (DESY, legacy):
  - CONF0: No Tungsten.
  - CONF[1,2,3]: Tungsten for different shower parts.
- [2021-11 (DESY)](https://twiki.cern.ch/twiki/bin/view/CALICE/SiWDESY202111):
  - CONF0: No Tungsten.
  - CONF1: 12 * 2.1 mm + 3 * 4.2 mm Tungsten.
- [2022-03 (DESY)](https://twiki.cern.ch/twiki/bin/viewauth/CALICE/SiWDESY202203):
  - CONF[0,2]: No Tungsten, different slab (Si thickness) arrangement (for different dates: 23/03 and 01/04).
  - CONF[1,3]: Same Tungsten configuration, different slab arrangement.
  - (For slab/ASU type arrangements, check [Layer arrangement](https://twiki.cern.ch/twiki/bin/viewauth/CALICE/SiWDESY202203#Layer_Arrangement).)
- 2022-06 (CERN):
  - CONF1: 7 x 2.8 mm + 8 x 4.2 mm Tungsten.
  - CONF2: 7 x 4.2 mm + 8 x 5.6 mm Tungsten.
  - CONF3: 24 x 4.2 mm Tungsten (ILD baseline 1). All 24 layers with 500 µm Si. 
  - CONF4: No Tungsten. Slab arrangement first proposed but not used in actual TB.
  - CONF5: 8 x 4.2 mm + 7 x 5.6 mm Tungsten, same slab arrangement as CONF4.
  - CONF6: 8 x 4.2 mm + 7 x 5.6 mm Tungsten. Slab arrangement used in june TB.
  - CONF7: No Tungsten, same slab arrangement as CONF6.
  - CONF[8,9]: Analogous to CONF[6,7] but handwritten to have 11mm squared cells for resolution studies.

### Visualization

Run 
```bash
geoDisplay generation/geometry/TB2022-03/ECAL_CONF0.xml # or any other geometry
```
There are some issues running this from LLR when one has zsh as default shell, but a workaround using x2go is:
```bash
source /cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/v02-00-01/init_ilcsoft.sh # still in zsh
bash
source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh # already in bash
```
Then run the `geoDisplay` as indicated above.

## Launching the simulation

The main script for preparing the simulation scripts (g4 macros, steering parameters, runnign shellscripts) and submitting the job to the condor queue is `send_condor.py`. The input of that script happens via `condor_config.yml` and the flags passed (explained in the help menu, `./send_condor.py --help`). A dry run (preparing without submitting) is possible by using `--run_locally`, then one could manually run the generation interactively or write a custom submit command.


## Comments

A few minor issues that need to be fixed:
- Some problems with the software configuration with zsh. Better to use bash.
- The script to write configurations is written in python3 (but ilcsoft gives python2). Needs to be run in a separate environment.

## Previous work

This part of the repository is based on [Daniel Jeans' work](https://gitlab.cern.ch/calice/calice_dd4heptestbeamsim/-/tree/master/2017_SiECAL_DESY/).
