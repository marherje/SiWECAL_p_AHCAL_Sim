## Note JESUS


TBD for ECAL+HCAL!






## Instruction to run the CED event display.
adrian.irles_at_ific.uv.es
2024/07/31

# DISCLAIMER

For some reason, the libglut3 libraries are not installed in glui01, but you cand download them from (last version) https://freeglut.sourceforge.net/docs/install.php and copy it to your local folder in glui01.ific.uv.es


# INSTALLATION

In your /lhome/ or wherever you decide:

> source /cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/init_ilcsoft.sh
> 
> tar xzvf freeglut-X.Y.Z.tar.gz (with X.Y.Z being the version of your download)
> 
> cd freeglut-X.Y.Z
> 
> mkdir build
> 
> cd build
> 
> cmake ..
> 
> make -j3
> 
> export LD_LIBRARY_PATH=${PWD}/freeglut-3.6.0/build/lib:$LD_LIBRARY_PATH
> 

--> this last step should be done everytime at the same time than source /cvmfs/ilc.desy.de/sw/x86_64_gcc103_
centos7/v02-03-03/init_ilcsoft.sh -> my  advice is that you make a local copy of        the init_ilcsoft.sh in your  own folder and include that line in the init_ilcsoft.sh   (you can run it every time you open a terminal by including it in your ~/.bashrc or just do it manually).

TAKE CARE! ${PWD} may need to be replaced by the full path !


# to make it run:

> cd Simplified_ECAL_PID_adrian/generation/run_scripts/ECALe_tests/test_display
> 
> glced &
> 
> Marlin display.xml
> 

--> NOTE the ECAL.xml used as argument in display.xml is a simple copy with of the ../../../geometry/ECALe_luxe/ECAL_CONF6_full.xml file with an extra line	 at the beginning:

> <global detectorName="ECALe_LUXE" />
> 
