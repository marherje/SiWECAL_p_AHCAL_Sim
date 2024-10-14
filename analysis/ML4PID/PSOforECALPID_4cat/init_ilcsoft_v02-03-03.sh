export ILCSOFT=/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03

# -------------------------------------------------------------------- ---

# ---  Use the same compiler and python as used for the installation   ---

# -------------------------------------------------------------------- ---
. /cvmfs/sft.cern.ch/lcg/releases/gcc/10.3.0-f5826/x86_64-centos7/setup.sh

export PATH=/cvmfs/sft.cern.ch/lcg/releases/Python/3.9.6-b0f98/x86_64-centos7-gcc10-opt/bin:${PATH}
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/Python/3.9.6-b0f98/x86_64-centos7-gcc10-opt/lib:${LD_LIBRARY_PATH}

export PYTHONPATH=/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21/python/examples:/cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc10-opt/lib/python3.9/site-packages

export CXX=g++
export CC=gcc
export GIT_EXEC_PATH=/cvmfs/sft.cern.ch/lcg/releases/git/2.29.2-e475b/x86_64-centos7-gcc10-opt/libexec/git-core

#--------------------------------------------------------------------------------
#     LCCD
#--------------------------------------------------------------------------------
export LCCD="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lccd/v01-05-02"
# --- additional LCCD commands ------- 


#--------------------------------------------------------------------------------
#     CondDBMySQL
#--------------------------------------------------------------------------------
export CondDBMySQL="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/CondDBMySQL/CondDBMySQL_ILC-0-9-7"
export COND_DB_DEBUGLOG="/dev/stdout"
# --- additional CondDBMySQL commands ------- 
export LD_LIBRARY_PATH="$CondDBMySQL/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCIO
#--------------------------------------------------------------------------------
export LCIO="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcio/v02-21"
# --- additional LCIO commands ------- 
export PATH="$LCIO/tools:$LCIO/bin:$PATH"
export LD_LIBRARY_PATH="$LCIO/lib:$LCIO/lib64:$LD_LIBRARY_PATH"
export PYTHONPATH="$LCIO/python:$LCIO/python/examples:$PYTHONPATH"
# --- additional ROOT commands ------- 
test -r /cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/root/6.30.04/bin/thisroot.sh && . /cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/root/6.30.04/bin/thisroot.sh


#--------------------------------------------------------------------------------
#     CMake
#--------------------------------------------------------------------------------
# --- additional CMake commands ------- 
export PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/CMake/3.28.3/bin:$PATH"


#--------------------------------------------------------------------------------
#     ILCUTIL
#--------------------------------------------------------------------------------
export ilcutil="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/ilcutil/v01-07-02"
# --- additional ILCUTIL commands ------- 
export LD_LIBRARY_PATH="$ilcutil/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Marlin
#--------------------------------------------------------------------------------
export MARLIN="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/Marlin/v01-19-02"
# --- additional Marlin commands ------- 
export PATH="$MARLIN/bin:$PATH"
export MARLIN_DLL="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/LCFIVertex/v00-09/lib/libLCFIVertexProcessors.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/MarlinDD4hep/v00-06-02/lib/libMarlinDD4hep.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/DDMarlinPandora/v00-12-01/lib/libDDMarlinPandora.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/MarlinKinfit/v00-06-01/lib/libMarlinKinfit.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/MarlinReco/v01-35/lib/libMarlinReco.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/PandoraAnalysis/v02-00-01/lib/libPandoraAnalysis.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/CEDViewer/v01-19-01/lib/libCEDViewer.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/Overlay/v00-23/lib/libOverlay.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/MarlinFastJet/v00-05-03/lib/libMarlinFastJet.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/LCTuple/v01-14/lib/libLCTuple.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/MarlinTrkProcessors/v02-12-05/lib/libMarlinTrkProcessors.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/MarlinKinfitProcessors/v00-05/lib/libMarlinKinfitProcessors.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/ILDPerformance/v01-12/lib/libILDPerformance.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/Clupatra/v01-03-01/lib/libClupatra.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/FCalClusterer/v01-00-06/lib/libFCalClusterer.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/LCFIPlus/v00-10-01/lib/libLCFIPlus.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/ForwardTracking/v01-14-02/lib/libForwardTracking.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/ConformalTracking/v01-12/lib/libConformalTracking.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/LICH/v00-01/lib/libLICH.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/Garlic/v03-01/lib/libGarlic.so:$MARLIN_DLL"


#--------------------------------------------------------------------------------
#     CLHEP
#--------------------------------------------------------------------------------
export CLHEP="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/CLHEP/2.4.7.1"
export CLHEP_BASE_DIR="$CLHEP"
export CLHEP_INCLUDE_DIR="$CLHEP/include"
# --- additional CLHEP commands ------- 
export PATH="$CLHEP_BASE_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$CLHEP_BASE_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     RAIDA
#--------------------------------------------------------------------------------
export RAIDA_HOME="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/RAIDA/v01-11"
# --- additional RAIDA commands ------- 
export PATH="$RAIDA_HOME/bin:$PATH"


#--------------------------------------------------------------------------------
#     GEAR
#--------------------------------------------------------------------------------
export GEAR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/gear/v01-09-02"
# --- additional GEAR commands ------- 
export PATH="$GEAR/tools:$GEAR/bin:$PATH"
export LD_LIBRARY_PATH="$GEAR/lib:$LD_LIBRARY_PATH"
# --- additional MarlinDD4hep commands ------- 


#--------------------------------------------------------------------------------
#     DD4hep
#--------------------------------------------------------------------------------
export DD4HEP="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/DD4hep/v01-28"
export DD4HEP_ENVINIT="${DD4HEP}/bin/thisdd4hep.sh"
export DD4hepExamples="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/DD4hep/v01-28/examples"
# --- additional DD4hep commands ------- 
test -r ${DD4HEP_ENVINIT} && . ${DD4HEP_ENVINIT}


#--------------------------------------------------------------------------------
#     Geant4
#--------------------------------------------------------------------------------
export G4INSTALL="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/geant4/11.2.1"
export G4SYSTEM="unknown-g++"
export G4ENV_INIT="$G4INSTALL/bin/geant4.sh"
# --- additional Geant4 commands ------- 
test -r ${G4ENV_INIT} && { cd $(dirname ${G4ENV_INIT}) ; . ./$(basename ${G4ENV_INIT}) ; cd $OLDPWD ; }


#--------------------------------------------------------------------------------
#     Qt5
#--------------------------------------------------------------------------------
export QTDIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/Qt5/v5.13.1"
# --- additional Qt5 commands ------- 
export PATH="$QTDIR/bin:$PATH"
export LD_LIBRARY_PATH="$QTDIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     XercesC
#--------------------------------------------------------------------------------
export XercesC_HOME="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/xercesc/v3.2.3"
# --- additional XercesC commands ------- 
export PATH="$XercesC_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$XercesC_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Boost
#--------------------------------------------------------------------------------
export BOOST_ROOT="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/boost/1.84.0"
# --- additional Boost commands ------- 


#--------------------------------------------------------------------------------
#     edm4hep
#--------------------------------------------------------------------------------
export EDM4HEP_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/edm4hep/v00-10-05/install"
# --- additional edm4hep commands ------- 
export LD_LIBRARY_PATH="$EDM4HEP_DIR/lib:$EDM4HEP_DIR/lib64:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     podio
#--------------------------------------------------------------------------------
export podio_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/podio/v00-99/install"
# --- additional podio commands ------- 
export PATH="$podio_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$podio_DIR/lib:$podio_DIR/lib64:$LD_LIBRARY_PATH"
export PYTHONPATH="$podio_DIR/python:$PYTHONPATH"
export ROOT_INCLUDE_PATH="$podio_DIR/include:$ROOT_INCLUDE_PATH"


#--------------------------------------------------------------------------------
#     SIO
#--------------------------------------------------------------------------------
export SIO_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/sio/v00-02"
# --- additional SIO commands ------- 
export PATH="$SIO_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$SIO_DIR/lib:$SIO_DIR/lib64:$LD_LIBRARY_PATH"
# --- additional DDKalTest commands ------- 


#--------------------------------------------------------------------------------
#     KalTest
#--------------------------------------------------------------------------------
export KALTEST="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/KalTest/v02-05-02"
# --- additional KalTest commands ------- 
export LD_LIBRARY_PATH="$KALTEST/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     aidaTT
#--------------------------------------------------------------------------------
export AIDATT="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/aidaTT/v00-10"
# --- additional aidaTT commands ------- 
export PATH="$AIDATT/bin:$PATH"
export LD_LIBRARY_PATH="$AIDATT/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     GSL
#--------------------------------------------------------------------------------
export GSL_HOME="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/gsl/2.7"
# --- additional GSL commands ------- 
export PATH="$GSL_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     GBL
#--------------------------------------------------------------------------------
export GBL="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/GBL/V02-02-01"
# --- additional GBL commands ------- 
export LD_LIBRARY_PATH="$GBL/lib:$LD_LIBRARY_PATH"
# --- additional Eigen commands ------- 
# --- additional DDMarlinPandora commands ------- 


#--------------------------------------------------------------------------------
#     MarlinUtil
#--------------------------------------------------------------------------------
# --- additional MarlinUtil commands ------- 
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/MarlinUtil/v01-17-02/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CED
#--------------------------------------------------------------------------------
# --- additional CED commands ------- 
export PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/CED/v01-09-04/bin:$PATH"
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/CED/v01-09-04/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     PandoraPFANew
#--------------------------------------------------------------------------------
export PANDORAPFANEW="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/PandoraPFANew/v03-25-03"
# --- additional PandoraPFANew commands ------- 
export LD_LIBRARY_PATH="$PANDORAPFANEW/lib:$LD_LIBRARY_PATH"
# --- additional MarlinTrk commands ------- 
# --- additional KalDet commands ------- 
# --- additional MarlinReco commands ------- 
# --- additional MarlinKinfit commands ------- 


#--------------------------------------------------------------------------------
#     PandoraAnalysis
#--------------------------------------------------------------------------------
export PANDORA_ANALYSIS_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/PandoraAnalysis/v02-00-01"
# --- additional PandoraAnalysis commands ------- 
export PATH="$PANDORA_ANALYSIS_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$PANDORA_ANALYSIS_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCFIVertex
#--------------------------------------------------------------------------------
# --- additional LCFIVertex commands ------- 
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/LCFIVertex/v00-09/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CEDViewer
#--------------------------------------------------------------------------------
# --- additional CEDViewer commands ------- 
export PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/CEDViewer/v01-19-01/bin:$PATH"
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/CEDViewer/v01-19-01/lib:$LD_LIBRARY_PATH"
# --- additional Overlay commands ------- 
# --- additional MarlinFastJet commands ------- 


#--------------------------------------------------------------------------------
#     FastJet
#--------------------------------------------------------------------------------
export FastJet_HOME="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/FastJet/3.4.2"
# --- additional FastJet commands ------- 
export PATH="$FastJet_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$FastJet_HOME/lib:$LD_LIBRARY_PATH"
# --- additional LCTuple commands ------- 
# --- additional KiTrack commands ------- 
# --- additional KiTrackMarlin commands ------- 
# --- additional MarlinTrkProcessors commands ------- 
# --- additional MarlinKinfitProcessors commands ------- 
# --- additional ILDPerformance commands ------- 
# --- additional Clupatra commands ------- 
# --- additional Physsim commands ------- 
# --- additional FCalClusterer commands ------- 
# --- additional LCFIPlus commands ------- 
# --- additional ForwardTracking commands ------- 
# --- additional ConformalTracking commands ------- 
# --- additional LICH commands ------- 
# --- additional Garlic commands ------- 


#--------------------------------------------------------------------------------
#     lcgeo
#--------------------------------------------------------------------------------
export lcgeo_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcgeo/v00-20-00"
export lcgeo_ENVINIT="${lcgeo_DIR}/bin/thislcgeo.sh"
export k4geo_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/lcgeo/v00-20-00"
# --- additional lcgeo commands ------- 
test -r ${lcgeo_ENVINIT} && . ${lcgeo_ENVINIT}
test -r ${k4geo_DIR}/bin/thisk4geo.sh && . ${k4geo_DIR}/bin/thisk4geo.sh


#--------------------------------------------------------------------------------
#     k4edm4hep2lcioconv
#--------------------------------------------------------------------------------
# --- additional k4edm4hep2lcioconv commands ------- 
export PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/k4edm4hep2lcioconv/v00-08-02/install/bin:$PATH"
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/k4edm4hep2lcioconv/v00-08-02/install/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/k4edm4hep2lcioconv/v00-08-02/install/lib64:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     MySQL
#--------------------------------------------------------------------------------
export MYSQL_HOME="/cvmfs/sft.cern.ch/lcg/releases/mysql/10.4.20-c0154/x86_64-centos7-gcc10-opt"
export MYSQL="$MYSQL_HOME"
export MYSQL_PATH="$MYSQL_HOME"
export MYSQL_LIBDIR="$MYSQL_HOME/lib64/mysql"
# --- additional MySQL commands ------- 
export PATH="$MYSQL_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$MYSQL_HOME/lib64:$MYSQL_HOME/lib:$MYSQL_HOME/lib64/mysql:$MYSQL_HOME/lib/mysql:$LD_LIBRARY_PATH"
