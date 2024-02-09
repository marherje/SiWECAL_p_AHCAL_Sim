#!/usr/bin/env python3

import argparse, sys, os, yaml
import subprocess, shlex
import numpy as np

# decide what's to be config-key and what's to be arg

parser = argparse.ArgumentParser(description="Send jobs to Condor queue.")
parser.add_argument('--particle', default="e-", help="Particle name, default \"e-\"")
parser.add_argument('--energy', default="3.0", type=str, help="Beam particle energy in GeV (default 3.0).")
parser.add_argument('--beam_x', default="0", type=str, help="Beam x position (deafult 0)")
parser.add_argument('--beam_y', default="0", type=str, help="Beam y position (deafult 0)")
parser.add_argument('--beam_sigma_x', default="7", type=str, help="Beam sigma x (deafult 7)")
parser.add_argument('--beam_sigma_y', default="7", type=str, help="Beam sigma y (deafult 7)")
parser.add_argument('--angle', default=0, help="Incidence angle of beam (testing this feature).")
parser.add_argument('--nevents', default=100, type=int, help="Number of events (default 100).")
parser.add_argument('--events_per_job', default=1000, type=int, help="Max number of events per batch job (default 1000).")
parser.add_argument('--first_job_index', default=1, type=int, help="Index of the first job/filename (default 1).")
parser.add_argument('--tbconf', default='TB2022-03_CONF0', help="TB+Conf as defined in confs.py (default TB2022-03_CONF0).")
# parser.add_argument('--test_run', action="store_true", default=False, help="Send test run to condor.")
parser.add_argument('--run_locally', action="store_true", help="Run locally and not in Condor.")

class Batch:
    def __init__(self,
                 particle,
                 energy,
                 beam_x,
                 beam_y,
                 beam_sigma_x,
                 beam_sigma_y,
                 angle,
                 nevents,
                 events_per_job,
                 first_job_index,
                 physics_list,
                 tbconf,
                 data_path,
                 ilcsoft_path,
                 geometry_folder,
                 submit_command,
                 run_locally): # run_locally is not used in the class

        self.particle = particle
        self.energy = energy
        self.beam_x = beam_x
        self.beam_y = beam_y
        self.beam_sigma_x = beam_sigma_x
        self.beam_sigma_y = beam_sigma_y
        self.angle = angle
        self.nevents = nevents 
        self.events_per_job = events_per_job 
        if nevents < events_per_job: self.events_per_job = nevents
        self.physics_list = physics_list
        self.tbconf = tbconf
        self.tb, self.conf = tbconf.split("_") 
        self.g4macroname = "{}/macros/{}_{}GeV.mac".format(self.tb, particle, energy)
        self.data_path = "{}{}/{}/lcio/".format(data_path, self.tb, self.conf)
        print("Data path is", data_path)
        self.cwd = os.getcwd()
        self.steer_path = "{}/{}/steer/".format(self.cwd, self.tb)
        self.geometry_folder = geometry_folder
        self.ilcsoft_path = ilcsoft_path
        self.submit_command = submit_command
        self.jobs_events = np.full(int(nevents / events_per_job)+1, events_per_job)
        if (nevents % events_per_job != 0):
            # rounding issue here can give 1 less event but...
            self.jobs_events[-1] = int(events_per_job * ((nevents / events_per_job) % 1))
        else:
            self.jobs_events = np.delete(self.jobs_events, -1)
        print("Launching", len(self.jobs_events), "jobs.\n", self.jobs_events)
        self.first_job_index = first_job_index
        # TODO: 
        #- Include mkdir CONF{N}/{build,lcio}
        # init counter of batch jobs


    def mkdirs(self):
        os.makedirs(self.data_path[:-5] + "build", exist_ok=True)
        os.makedirs(self.data_path, exist_ok=True)

    def write_g4macro(self):
        for it, nev in enumerate(self.jobs_events): 
            ijob = it + self.first_job_index 
            macroname = self.g4macroname[:-4] + "_" + str(ijob) + ".mac"
            with open(macroname, 'w') as f:
                f.write("/gps/verbose 1\n")
                f.write("/gps/particle {}\n".format(self.particle))
                f.write("/gps/direction 0 0 1\n")
                f.write("/gps/pos/type Beam\n")
                f.write("/gps/pos/shape Circle\n")
                f.write("/gps/pos/centre {} {} 0 mm\n".format(self.beam_x, self.beam_y))
                f.write("/gps/pos/sigma_x {} mm\n".format(self.beam_sigma_x))
                f.write("/gps/pos/sigma_y {} mm\n".format(self.beam_sigma_y))
                # Angle should be added here
                f.write("/gps/ang/rot1 0 0 1\n")
                f.write("/gps/ang/rot2 0 1 0\n")
                f.write("/gps/ene/type Mono\n")
                f.write("/gps/ene/mono {} GeV\n".format(self.energy))
                f.write("/run/beamOn {}".format(nev))
        print("G4 macro written in", (self.g4macroname))

    
    def write_pyscripts(self):
        for it, _ in enumerate(self.jobs_events):
            ijob = it + self.first_job_index
            label = "{}_{}_{}_{}GeV_{}".format(self.physics_list,
             self.tbconf,
             self.particle,
             self.energy,
             ijob)  
             # Should be only one label...
            label2 = "ECAL_{}_conf{}_{}_{}GeV_{}".format(self.physics_list,
                                                         #self.conf[-1], ## Sould be debugged!!!
                                                         self.conf.split("CONF")[-1], ## testing!!!
                                                         self.particle,
                                                         self.energy,
                                                         ijob)
            print("Label is", label)
            pyscript = self.steer_path + "/runddsim_" + label + ".py"
            with open(pyscript, 'w') as f:
                f.write("from DDSim.DD4hepSimulation import  DD4hepSimulation\n")
                f.write("#from SystemOfUnits import mm, GeV, MeV\n")
                f.write("from g4units import GeV, mm, MeV\n")
                f.write("\n")
                f.write("SIM = DD4hepSimulation()\n")
                f.write("\n")
                f.write("SIM.runType = \"run\"\n")
                f.write("# Number of events defined in macro file\n")
                f.write("\n")
                f.write("SIM.skipNEvents = 0\n")
                f.write("SIM.outputFile = \"{}{}.slcio\"\n".format(self.data_path, label2))
                f.write("\n")
                f.write("SIM.compactFile = \"{}/{}/ECAL_{}.xml\"\n".format(self.geometry_folder, self.tb, self.conf))
                f.write("SIM.dumpSteeringFile = \"{}dumpSteering.xml\"\n".format(self.steer_path))
                f.write("\n")
                f.write("SIM.field.eps_min = 1*mm\n")
                f.write("SIM.part.minimalKineticEnergy = 0.3*MeV\n")
                f.write("SIM.physicsList = \"{}\"\n".format(self.physics_list))
                f.write("SIM.enableDetailedShowerMode=True")
            print("Python script written in", pyscript)
            # break
    
    def write_shcondors(self):
        #for it in range(1, 21):
        for it, _ in enumerate(self.jobs_events):
            ijob = it + self.first_job_index
            label = "{}_{}_{}_{}GeV_{}".format(self.physics_list,
             self.tbconf,
             self.particle,
             self.energy,
             ijob)  
            pyscript = self.steer_path + "runddsim_" + label + ".py"
            condorsh = self.steer_path + "/runddsim_" + label + ".sh"
            macroname = self.g4macroname[:-4] + "_" + str(ijob) + ".mac"
            with open(condorsh, 'w') as f:
                f.write("#!/bin/sh\n")
                f.write("source {}/init_ilcsoft.sh\n".format(self.ilcsoft_path)) 
                print(os.getcwd())
                f.write("cp -r {}/runddsim_{}.{{py,sh}} .\n".format(self.steer_path, label))
                f.write("#This is run in {}\n".format(os.getcwd()))
                #f.write("ddsim --enableG4GPS --macroFile {}/{} --steeringFile {}\n".format(os.getcwd(), self.g4macroname, pyscript))
                f.write("ddsim --enableG4GPS --macroFile {}/{} --steeringFile {}\n".format(os.getcwd(), macroname, pyscript))
            os.chmod(condorsh, 0o755)
            print("Condor script written in", condorsh)

    def submit_jobs(self):
        for it, _ in enumerate(self.jobs_events):
            ijob = it + self.first_job_index
            label = "{}_{}_{}_{}GeV_{}".format(self.physics_list,
             self.tbconf,
             self.particle,
             self.energy,
             ijob)
            #condorsh = self.steer_path + "/runddsim_" + label + ".sh"
            condorsh = "runddsim_" + label + ".sh"
            wd = os.getcwd()
            os.chdir(self.steer_path)
            #command = "/opt/exp_soft/cms/t3/t3submit -short {}".format(condorsh)
            command = self.submit_command + " " + condorsh
            shargs = shlex.split(command)
            p = subprocess.Popen(shargs)
            os.chdir(wd)
            print(command)
            # break

if __name__ == "__main__":
    args = parser.parse_args()
    with open("condor_config.yml", "r") as ymlfile:
        config = yaml.safe_load(ymlfile)
    
    batch = Batch(**vars(args), **config)
    print(batch.steer_path)
    batch.mkdirs()
    batch.write_g4macro()
    batch.write_pyscripts()
    batch.write_shcondors()
    if not args.run_locally: batch.submit_jobs()
    else: print("Run locally condor script!!")
