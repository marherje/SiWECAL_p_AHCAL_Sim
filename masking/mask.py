import os
import sys
import glob
import subprocess
import argparse
from pathlib import Path
currentDir = Path(__file__).parent.absolute()

simDir = os.path.join(currentDir, "Data", "sim")
RunSettingsDir = os.path.join(currentDir, "Run_Settings")

runIDs = {
  10 : 90320,
  20 : 90378,
  40 : 90375,
  60 : 90372,
  80 : 90367,
  100 : 90365,
  150 : 90355
}

for ienergy in [10, 20, 40, 60, 80, 100, 150]:

  simFile = os.path.join(simDir, f'ECAL_QGSP_BERT_conf6_e-_{ienergy}.0GeV_build.root')

  runID = runIDs[ienergy]
  RunSettings = os.path.join(RunSettingsDir, f'Run_Settings_{runID}_e-_{ienergy}.0GeV.txt')

  cmd = f'root -l masker.C\\(\\\"{RunSettings}\\\",\\\"{simFile}\\\"\\)'
  print(cmd)
  subprocess.run(cmd, shell=True)