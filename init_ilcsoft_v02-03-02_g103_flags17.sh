export ILCSOFT=/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02

# -------------------------------------------------------------------- ---

# ---  Use the same compiler and python as used for the installation   ---

# -------------------------------------------------------------------- ---
. /cvmfs/sft.cern.ch/lcg/releases/gcc/10.3.0-f5826/x86_64-centos7/setup.sh

export PATH=/cvmfs/sft.cern.ch/lcg/releases/Python/3.9.6-b0f98/x86_64-centos7-gcc10-opt/bin:${PATH}
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/Python/3.9.6-b0f98/x86_64-centos7-gcc10-opt/lib:${LD_LIBRARY_PATH}

export PYTHONPATH=/cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc10-opt/lib/python3.9/site-packages

export CXX=g++
export CC=gcc
export GIT_EXEC_PATH=/cvmfs/sft.cern.ch/lcg/releases/git/2.29.2-e475b/x86_64-centos7-gcc10-opt/libexec/git-core
export CXXFLAGS="-std=c++17"

#--------------------------------------------------------------------------------
#     LCCD
#--------------------------------------------------------------------------------
export LCCD="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/lccd/v01-05-01"
# --- additional LCCD commands ------- 


#--------------------------------------------------------------------------------
#     CondDBMySQL
#--------------------------------------------------------------------------------
export CondDBMySQL="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/CondDBMySQL/CondDBMySQL_ILC-0-9-7"
export COND_DB_DEBUGLOG="/dev/stdout"
# --- additional CondDBMySQL commands ------- 
export LD_LIBRARY_PATH="$CondDBMySQL/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCIO
#--------------------------------------------------------------------------------
export LCIO="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/lcio/v02-20"
# --- additional LCIO commands ------- 
export PATH="$LCIO/tools:$LCIO/bin:$PATH"
export LD_LIBRARY_PATH="$LCIO/lib:$LCIO/lib64:$LD_LIBRARY_PATH"
export PYTHONPATH="$LCIO/python:$LCIO/python/examples:$PYTHONPATH"
# --- additional ROOT commands ------- 
test -r /cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/root/6.28.04/bin/thisroot.sh && . /cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/root/6.28.04/bin/thisroot.sh


#--------------------------------------------------------------------------------
#     Marlin
#--------------------------------------------------------------------------------
export MARLIN="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/Marlin/v01-19"
# --- additional Marlin commands ------- 
export PATH="$MARLIN/bin:$PATH"
export MARLIN_DLL="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/MarlinDD4hep/v00-06-02/lib/libMarlinDD4hep.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/DDMarlinPandora/v00-12/lib/libDDMarlinPandora.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/MarlinReco/v01-34/lib/libMarlinReco.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/PandoraAnalysis/v02-00-01/lib/libPandoraAnalysis.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/LCFIVertex/v00-08/lib/libLCFIVertexProcessors.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/CEDViewer/v01-19-01/lib/libCEDViewer.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/Overlay/v00-23/lib/libOverlay.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/MarlinFastJet/v00-05-03/lib/libMarlinFastJet.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/LCTuple/v01-14/lib/libLCTuple.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/MarlinKinfit/v00-06-01/lib/libMarlinKinfit.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/MarlinTrkProcessors/v02-12-03/lib/libMarlinTrkProcessors.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/MarlinKinfitProcessors/v00-05/lib/libMarlinKinfitProcessors.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/ILDPerformance/v01-12/lib/libILDPerformance.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/Clupatra/v01-03/lib/libClupatra.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/Physsim/v00-04-02/lib/libPhyssim.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/FCalClusterer/v01-00-03/lib/libFCalClusterer.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/LCFIPlus/v00-10-01/lib/libLCFIPlus.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/ForwardTracking/v01-14-01/lib/libForwardTracking.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/ConformalTracking/v01-11-01/lib/libConformalTracking.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/LICH/v00-01/lib/libLICH.so:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/Garlic/v03-01/lib/libGarlic.so:$MARLIN_DLL"


#--------------------------------------------------------------------------------
#     CLHEP
#--------------------------------------------------------------------------------
export CLHEP="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/CLHEP/2.4.6.4"
export CLHEP_BASE_DIR="$CLHEP"
export CLHEP_INCLUDE_DIR="$CLHEP/include"
# --- additional CLHEP commands ------- 
export PATH="$CLHEP_BASE_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$CLHEP_BASE_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     RAIDA
#--------------------------------------------------------------------------------
export RAIDA_HOME="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/RAIDA/v01-11"
# --- additional RAIDA commands ------- 
export PATH="$RAIDA_HOME/bin:$PATH"


#--------------------------------------------------------------------------------
#     GEAR
#--------------------------------------------------------------------------------
export GEAR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/gear/v01-09-01"
# --- additional GEAR commands ------- 
export PATH="$GEAR/tools:$GEAR/bin:$PATH"
export LD_LIBRARY_PATH="$GEAR/lib:$LD_LIBRARY_PATH"
# --- additional MarlinDD4hep commands ------- 


#--------------------------------------------------------------------------------
#     DD4hep
#--------------------------------------------------------------------------------
export DD4HEP="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/DD4hep/v01-25-01"
export DD4HEP_ENVINIT="${DD4HEP}/bin/thisdd4hep.sh"
# --- additional DD4hep commands ------- 
test -r ${DD4HEP_ENVINIT} && . ${DD4HEP_ENVINIT}


#--------------------------------------------------------------------------------
#     Geant4
#--------------------------------------------------------------------------------
export G4INSTALL="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/geant4/11.1.1"
export G4SYSTEM="unknown-g++"
export G4ENV_INIT="$G4INSTALL/bin/geant4.sh"
# --- additional Geant4 commands ------- 
test -r ${G4ENV_INIT} && { cd $(dirname ${G4ENV_INIT}) ; . ./$(basename ${G4ENV_INIT}) ; cd $OLDPWD ; }


#--------------------------------------------------------------------------------
#     Qt5
#--------------------------------------------------------------------------------
export QTDIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/Qt5/v5.13.1"
# --- additional Qt5 commands ------- 
export PATH="$QTDIR/bin:$PATH"
export LD_LIBRARY_PATH="$QTDIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     XercesC
#--------------------------------------------------------------------------------
export XercesC_HOME="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/xercesc/v3.2.3"
# --- additional XercesC commands ------- 
export PATH="$XercesC_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$XercesC_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Boost
#--------------------------------------------------------------------------------
export BOOST_ROOT="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/boost/1.77.0"
# --- additional Boost commands ------- 


#--------------------------------------------------------------------------------
#     edm4hep
#--------------------------------------------------------------------------------
export EDM4HEP_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/edm4hep/v00-10/install"
# --- additional edm4hep commands ------- 
export LD_LIBRARY_PATH="$EDM4HEP_DIR/lib:$EDM4HEP_DIR/lib64:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     podio
#--------------------------------------------------------------------------------
export podio_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/podio/v00-16-06/install"
# --- additional podio commands ------- 
export PATH="$podio_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$podio_DIR/lib:$podio_DIR/lib64:$LD_LIBRARY_PATH"
export PYTHONPATH="$podio_DIR/python:$PYTHONPATH"
export ROOT_INCLUDE_PATH="$podio_DIR/include:$ROOT_INCLUDE_PATH"


#--------------------------------------------------------------------------------
#     SIO
#--------------------------------------------------------------------------------
export SIO_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/sio/v00-01"
# --- additional SIO commands ------- 
export PATH="$SIO_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$SIO_DIR/lib:$SIO_DIR/lib64:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CMake
#--------------------------------------------------------------------------------
# --- additional CMake commands ------- 
export PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/CMake/3.23.2/bin:$PATH"
# --- additional DDKalTest commands ------- 


#--------------------------------------------------------------------------------
#     KalTest
#--------------------------------------------------------------------------------
export KALTEST="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/KalTest/v02-05-01"
# --- additional KalTest commands ------- 
export LD_LIBRARY_PATH="$KALTEST/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     aidaTT
#--------------------------------------------------------------------------------
export AIDATT="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/aidaTT/v00-10"
# --- additional aidaTT commands ------- 
export PATH="$AIDATT/bin:$PATH"
export LD_LIBRARY_PATH="$AIDATT/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     GSL
#--------------------------------------------------------------------------------
export GSL_HOME="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/gsl/2.7"
# --- additional GSL commands ------- 
export PATH="$GSL_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     GBL
#--------------------------------------------------------------------------------
export GBL="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/GBL/V02-02-01"
# --- additional GBL commands ------- 
export LD_LIBRARY_PATH="$GBL/lib:$LD_LIBRARY_PATH"
# --- additional Eigen commands ------- 
# --- additional DDMarlinPandora commands ------- 


#--------------------------------------------------------------------------------
#     MarlinUtil
#--------------------------------------------------------------------------------
# --- additional MarlinUtil commands ------- 
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/MarlinUtil/v01-17-01/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CED
#--------------------------------------------------------------------------------
# --- additional CED commands ------- 
export PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/CED/v01-09-04/bin:$PATH"
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/CED/v01-09-04/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     PandoraPFANew
#--------------------------------------------------------------------------------
export PANDORAPFANEW="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/PandoraPFANew/v03-25-03"
# --- additional PandoraPFANew commands ------- 
export LD_LIBRARY_PATH="$PANDORAPFANEW/lib:$LD_LIBRARY_PATH"
# --- additional MarlinTrk commands ------- 
# --- additional KalDet commands ------- 
# --- additional MarlinReco commands ------- 
# --- additional MarlinKinfit commands ------- 


#--------------------------------------------------------------------------------
#     PandoraAnalysis
#--------------------------------------------------------------------------------
export PANDORA_ANALYSIS_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/PandoraAnalysis/v02-00-01"
# --- additional PandoraAnalysis commands ------- 
export PATH="$PANDORA_ANALYSIS_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$PANDORA_ANALYSIS_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCFIVertex
#--------------------------------------------------------------------------------
# --- additional LCFIVertex commands ------- 
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/LCFIVertex/v00-08/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CEDViewer
#--------------------------------------------------------------------------------
# --- additional CEDViewer commands ------- 
export PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/CEDViewer/v01-19-01/bin:$PATH"
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/CEDViewer/v01-19-01/lib:$LD_LIBRARY_PATH"
# --- additional Overlay commands ------- 
# --- additional MarlinFastJet commands ------- 


#--------------------------------------------------------------------------------
#     FastJet
#--------------------------------------------------------------------------------
export FastJet_HOME="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/FastJet/3.4.0"
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
export lcgeo_DIR="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/lcgeo/v00-18-01"
export lcgeo_ENVINIT="${lcgeo_DIR}/bin/thislcgeo.sh"
# --- additional lcgeo commands ------- 
test -r ${lcgeo_ENVINIT} && . ${lcgeo_ENVINIT}


#--------------------------------------------------------------------------------
#     k4edm4hep2lcioconv
#--------------------------------------------------------------------------------
# --- additional k4edm4hep2lcioconv commands ------- 
export PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/k4edm4hep2lcioconv/v00-05/install/bin:$PATH"
export LD_LIBRARY_PATH="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/k4edm4hep2lcioconv/v00-05/install/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/k4edm4hep2lcioconv/v00-05/install/lib64:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     DD4hepExamples
#--------------------------------------------------------------------------------
export DD4hepExamples="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/DD4hepExamples/v01-25-01"
# --- additional DD4hepExamples commands ------- 
export PATH="$DD4hepExamples/bin:$PATH"
export LD_LIBRARY_PATH="$DD4hepExamples/lib:$LD_LIBRARY_PATH"


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


#--------------------------------------------------------------------------------
#     ILCUTIL
#--------------------------------------------------------------------------------
export ilcutil="/cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-02/ilcutil/v01-07-01"
# --- additional ILCUTIL commands ------- 
export LD_LIBRARY_PATH="$ilcutil/lib:$LD_LIBRARY_PATH"
