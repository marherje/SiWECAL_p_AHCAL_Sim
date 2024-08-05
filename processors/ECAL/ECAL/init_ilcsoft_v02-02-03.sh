export ILCSOFT=/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03

# -------------------------------------------------------------------- ---

# ---  Use the same compiler and python as used for the installation   ---

# -------------------------------------------------------------------- ---
. /cvmfs/sft.cern.ch/lcg/releases/gcc/8.2.0-3fa06/x86_64-centos7/setup.sh

export PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_96/Python/2.7.16/x86_64-centos7-gcc8-opt/bin:${PATH}
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_96/Python/2.7.16/x86_64-centos7-gcc8-opt/lib:${LD_LIBRARY_PATH}

export CXX=g++
export CC=gcc


#--------------------------------------------------------------------------------
#     LCCD
#--------------------------------------------------------------------------------
export LCCD="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/lccd/v01-05"
# --- additional LCCD commands ------- 


#--------------------------------------------------------------------------------
#     CondDBMySQL
#--------------------------------------------------------------------------------
export COND_DB_DEBUGLOG="/dev/stdout"
export CondDBMySQL="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CondDBMySQL/CondDBMySQL_ILC-0-9-7"
# --- additional CondDBMySQL commands ------- 
export LD_LIBRARY_PATH="$CondDBMySQL/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCIO
#--------------------------------------------------------------------------------
export LCIO="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/lcio/v02-17"
# --- additional LCIO commands ------- 
export PYTHONPATH="$LCIO/python:$LCIO/python/examples:$PYTHONPATH"
export PATH="$LCIO/tools:$LCIO/bin:$PATH"
export LD_LIBRARY_PATH="$LCIO/lib:$LCIO/lib64:$LD_LIBRARY_PATH"
# --- additional ROOT commands ------- 
test -r /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/bin/thisroot.sh && . /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/bin/thisroot.sh


#--------------------------------------------------------------------------------
#     CMake
#--------------------------------------------------------------------------------
# --- additional CMake commands ------- 
export PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/CMake/3.15.5/bin:$PATH"


#--------------------------------------------------------------------------------
#     ILCUTIL
#--------------------------------------------------------------------------------
export ilcutil="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/ilcutil/v01-06-02"
# --- additional ILCUTIL commands ------- 
export LD_LIBRARY_PATH="$ilcutil/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Marlin
#--------------------------------------------------------------------------------
export MARLIN="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Marlin/v01-17-01"
# --- additional Marlin commands ------- 
export PATH="$MARLIN/bin:$PATH"
export MARLIN_DLL="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinDD4hep/v00-06/lib/libMarlinDD4hep.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/DDMarlinPandora/v00-12/lib/libDDMarlinPandora.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinReco/v01-32/lib/libMarlinReco.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/PandoraAnalysis/v02-00-01/lib/libPandoraAnalysis.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/LCFIVertex/v00-08/lib/libLCFIVertexProcessors.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CEDViewer/v01-19/lib/libCEDViewer.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Overlay/v00-22-03/lib/libOverlay.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinFastJet/v00-05-02/lib/libMarlinFastJet.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/LCTuple/v01-13/lib/libLCTuple.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinKinfit/v00-06/lib/libMarlinKinfit.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinTrkProcessors/v02-12/lib/libMarlinTrkProcessors.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinKinfitProcessors/v00-04-02/lib/libMarlinKinfitProcessors.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/ILDPerformance/v01-10/lib/libILDPerformance.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Clupatra/v01-03/lib/libClupatra.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Physsim/v00-04-01/lib/libPhyssim.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/LCFIPlus/v00-10-01/lib/libLCFIPlus.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/FCalClusterer/v01-00-01/lib/libFCalClusterer.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/ForwardTracking/v01-14/lib/libForwardTracking.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/ConformalTracking/v01-11/lib/libConformalTracking.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/LICH/v00-01/lib/libLICH.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Garlic/v03-01/lib/libGarlic.so:$MARLIN_DLL"


#--------------------------------------------------------------------------------
#     CLHEP
#--------------------------------------------------------------------------------
export CLHEP="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CLHEP/2.3.4.3"
export CLHEP_BASE_DIR="$CLHEP"
export CLHEP_INCLUDE_DIR="$CLHEP/include"
# --- additional CLHEP commands ------- 
export PATH="$CLHEP_BASE_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$CLHEP_BASE_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     RAIDA
#--------------------------------------------------------------------------------
export RAIDA_HOME="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/RAIDA/v01-09"
# --- additional RAIDA commands ------- 
export PATH="$RAIDA_HOME/bin:$PATH"


#--------------------------------------------------------------------------------
#     GEAR
#--------------------------------------------------------------------------------
export GEAR="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09"
# --- additional GEAR commands ------- 
export PATH="$GEAR/tools:$GEAR/bin:$PATH"
export LD_LIBRARY_PATH="$GEAR/lib:$LD_LIBRARY_PATH"
# --- additional MarlinDD4hep commands ------- 


#--------------------------------------------------------------------------------
#     DD4hep
#--------------------------------------------------------------------------------
export DD4HEP="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/DD4hep/v01-11-02"
export DD4HEP_ENVINIT="${DD4HEP}/bin/thisdd4hep.sh"
# --- additional DD4hep commands ------- 
test -r ${DD4HEP_ENVINIT} && . ${DD4HEP_ENVINIT}


#--------------------------------------------------------------------------------
#     Geant4
#--------------------------------------------------------------------------------
export G4INSTALL="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/geant4/10.04.p03"
export G4ENV_INIT="$G4INSTALL/bin/geant4.sh"
export G4SYSTEM="Linux-g++"
# --- additional Geant4 commands ------- 
test -r ${G4ENV_INIT} && { cd $(dirname ${G4ENV_INIT}) ; . ./$(basename ${G4ENV_INIT}) ; cd $OLDPWD ; }


#--------------------------------------------------------------------------------
#     Qt5
#--------------------------------------------------------------------------------
export QTDIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Qt5/v5.13.1"
# --- additional Qt5 commands ------- 
export PATH="$QTDIR/bin:$PATH"
export LD_LIBRARY_PATH="$QTDIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     XercesC
#--------------------------------------------------------------------------------
export XercesC_HOME="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/xercesc/v3.2.2"
# --- additional XercesC commands ------- 
export PATH="$XercesC_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$XercesC_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Boost
#--------------------------------------------------------------------------------
export BOOST_ROOT="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/boost/1.71.0"
# --- additional Boost commands ------- 
# --- additional DDKalTest commands ------- 


#--------------------------------------------------------------------------------
#     KalTest
#--------------------------------------------------------------------------------
export KALTEST="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/KalTest/v02-05"
# --- additional KalTest commands ------- 
export LD_LIBRARY_PATH="$KALTEST/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     aidaTT
#--------------------------------------------------------------------------------
export AIDATT="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/aidaTT/v00-10"
# --- additional aidaTT commands ------- 
export PATH="$AIDATT/bin:$PATH"
export LD_LIBRARY_PATH="$AIDATT/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     GSL
#--------------------------------------------------------------------------------
export GSL_HOME="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gsl/2.6"
# --- additional GSL commands ------- 
export PATH="$GSL_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     GBL
#--------------------------------------------------------------------------------
export GBL="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/GBL/V02-02-00"
# --- additional GBL commands ------- 
export LD_LIBRARY_PATH="$GBL/lib:$LD_LIBRARY_PATH"
# --- additional Eigen commands ------- 
# --- additional DDMarlinPandora commands ------- 


#--------------------------------------------------------------------------------
#     MarlinUtil
#--------------------------------------------------------------------------------
# --- additional MarlinUtil commands ------- 
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/MarlinUtil/v01-16-01/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CED
#--------------------------------------------------------------------------------
# --- additional CED commands ------- 
export PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CED/v01-09-03/bin:$PATH"
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CED/v01-09-03/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     PandoraPFANew
#--------------------------------------------------------------------------------
export PANDORAPFANEW="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/PandoraPFANew/v03-13-02"
# --- additional PandoraPFANew commands ------- 
export LD_LIBRARY_PATH="$PANDORAPFANEW/lib:$LD_LIBRARY_PATH"
# --- additional MarlinTrk commands ------- 
# --- additional KalDet commands ------- 
# --- additional MarlinReco commands ------- 
# --- additional MarlinKinfit commands ------- 


#--------------------------------------------------------------------------------
#     PandoraAnalysis
#--------------------------------------------------------------------------------
export PANDORA_ANALYSIS_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/PandoraAnalysis/v02-00-01"
# --- additional PandoraAnalysis commands ------- 
export PATH="$PANDORA_ANALYSIS_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$PANDORA_ANALYSIS_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCFIVertex
#--------------------------------------------------------------------------------
# --- additional LCFIVertex commands ------- 
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/LCFIVertex/v00-08/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CEDViewer
#--------------------------------------------------------------------------------
# --- additional CEDViewer commands ------- 
export PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CEDViewer/v01-19/bin:$PATH"
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CEDViewer/v01-19/lib:$LD_LIBRARY_PATH"
# --- additional Overlay commands ------- 
# --- additional MarlinFastJet commands ------- 


#--------------------------------------------------------------------------------
#     FastJet
#--------------------------------------------------------------------------------
export FastJet_HOME="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/FastJet/3.2.1"
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
# --- additional LCFIPlus commands ------- 
# --- additional FCalClusterer commands ------- 
# --- additional ForwardTracking commands ------- 
# --- additional ConformalTracking commands ------- 
# --- additional LICH commands ------- 
# --- additional PathFinder commands ------- 
# --- additional BBQ commands ------- 
# --- additional Garlic commands ------- 


#--------------------------------------------------------------------------------
#     lcgeo
#--------------------------------------------------------------------------------
export lcgeo_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/lcgeo/v00-16-07"
export lcgeo_ENVINIT="${lcgeo_DIR}/bin/thislcgeo.sh"
# --- additional lcgeo commands ------- 
test -r ${lcgeo_ENVINIT} && . ${lcgeo_ENVINIT}


#--------------------------------------------------------------------------------
#     DD4hepExamples
#--------------------------------------------------------------------------------
export DD4hepExamples="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/DD4hepExamples/v01-11-02"
# --- additional DD4hepExamples commands ------- 
export PATH="$DD4hepExamples/bin:$PATH"
export LD_LIBRARY_PATH="$DD4hepExamples/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     SIO
#--------------------------------------------------------------------------------
export SIO_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/sio/v00-01"
# --- additional SIO commands ------- 
export PATH="$SIO_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$SIO_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     MySQL
#--------------------------------------------------------------------------------
export MYSQL_HOME="/cvmfs/sft.cern.ch/lcg/releases/mysql/5.7.26-c3e26/x86_64-centos7-gcc8-opt"
export MYSQL_LIBDIR="$MYSQL_HOME/lib64/mysql"
export MYSQL_PATH="$MYSQL_HOME"
export MYSQL="$MYSQL_HOME"
# --- additional MySQL commands ------- 
export PATH="$MYSQL_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$MYSQL_HOME/lib64:$MYSQL_HOME/lib:$MYSQL_HOME/lib64/mysql:$MYSQL_HOME/lib/mysql:$LD_LIBRARY_PATH"
