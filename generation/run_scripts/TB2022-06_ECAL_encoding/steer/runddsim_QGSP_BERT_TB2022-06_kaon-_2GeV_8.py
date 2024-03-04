from DDSim.DD4hepSimulation import DD4hepSimulation
#from SystemOfUnits import mm, GeV, MeV
from g4units import GeV, mm, MeV

SIM = DD4hepSimulation()

SIM.runType = "run"
# Number of events defined in macro file
#SIM.numberOfEvents = 3000

SIM.skipNEvents = 0
SIM.outputFile = "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/generation/run_scripts/TB2022-06_ECAL_encoding/data/TB2022-06_QGSP_BERT_TB2022-06_kaon-_2GeV_8.slcio"
# SIM.outputFile = "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/generation/run_scripts/TB2022-06_ECAL_encoding/data/TB2022-06_QGSP_BERT_TB2022-06_kaon-_2GeV_8.slcio"

SIM.compactFile = "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/generation/geometry/ECALHCAL_TB2022-06/TBModel_June2022.xml"
SIM.dumpSteeringFile = "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/generation/run_scripts/TB2022-06_ECAL_encoding/steer/dumpSteering.xml"

SIM.field.eps_min = 0.0001*mm
SIM.part.minimalKineticEnergy = 0.3*MeV
SIM.physicsList = "QGSP_BERT"
SIM.enableDetailedShowerMode=True
