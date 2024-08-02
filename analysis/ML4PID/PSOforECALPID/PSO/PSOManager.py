#!/us/bin/env python
# -*- coding: utf-8 -*-
import os, sys, json, math, ROOT, time, subprocess

from particle  import Particle
from QueHelper import QueHelper

from common import *

class PSOManager:
    def __init__(self, OutputDir='', DataSubdir='InitData', Verbose=True, PSOConfig=''):
      self.Verbose=Verbose
      self.vw=1.0
      self.vp=1.0
      self.vg=1.0
      self.Particles=[]
      self.BestFOMGlobal=0.0
      self.BestADGlobal=0.0
      self.BestKSGlobal=0.0
      self.BestCoordinatesGlobal=[]
      self.Coordinates=[]

      self.rand=ROOT.TRandom3(0) ## change to !=0 if you want repeatable random seeds

      self.OutputDir  = os.path.abspath('.' if OutputDir  == '' else OutputDir)
      self.DataSubdir = self.OutputDir+'/'+('InitData' if DataSubdir == '' else DataSubdir)

      subprocess.call(['mkdir', '-p', self.DataSubdir+'/logs'])
      subprocess.call(['mkdir', '-p', self.DataSubdir+'/weights'])

      subprocess.call(['touch', self.DataSubdir+'/ParticleConfig.txt'])
      subprocess.call(['touch', self.DataSubdir+'/ParticleResult.txt'])
      subprocess.call(['touch', self.DataSubdir+'/autoLog.txt'])
      subprocess.call(['touch', self.DataSubdir+'/dump.txt'])

      self.nParticles=10
      self.TenBestMVAs=[[0.0,0.0,0.0,"",[],[],[]] for i in range(10)]
      #print self.TenBestMVAs
      self.ADThreshold=0.1
      self.KSThreshold=0.01
      
      RunSystem = 'IFIC'

      self.NIterations=1
      self.FOM=""
      self.FindBestVariables = 1
      self.MaxVariablesInCombination=10
      self.ImprovementThreshold=1.00
      self.FactoryString=""
      self.PreparationString=""
      self.BWeightExpression=""
      self.CWeightExpression=""
      self.OWeightExpression=""
      self.MethodType=""
      self.MethodParams=""
      self.BTreeName=""
      self.CTreeName=""
      self.OTreeName=""
      self.RepeatTrainingNTimes=0
      self.DrawNRandomAsStartingVars=0
      self.SaveTrainingsToTrees="False"
      self.UseFixedTrainTestSplitting = 0
      self.UseFixedTrainTestSplitting_Train = '(eventNumber%2)==0'

      initialVariables = []
      additionalVariables = []
      SpectatorVariables = []
      ivars = []
      avars = []
      svars = []
      weightVariables = []
      configFile=open(PSOConfig,"r")
      configLines=list(configFile)
      configFile.close()
      for i,l in enumerate(configLines):
        configLines[i]=configLines[i].replace("\n","")

      if "InitialVariables:" in configLines and "EndVariables" in configLines:
        id1=configLines.index("InitialVariables:")
        id2=configLines.index("EndVariables")
        initialVariables=configLines[id1+1:id2]
        configLines=configLines[:id1]+configLines[id2+1:]
      if "AdditionalVariables:" in configLines and "EndVariables" in configLines:
        id1=configLines.index("AdditionalVariables:")
        id2=configLines.index("EndVariables")
        additionalVariables=configLines[id1+1:id2]
        configLines=configLines[:id1]+configLines[id2+1:]
      if "SpectatorVariables:" in configLines and "EndVariables" in configLines:
        id1=configLines.index("SpectatorVariables:")
        id2=configLines.index("EndVariables")
        SpectatorVariables=configLines[id1+1:id2]
        configLines=configLines[:id1]+configLines[id2+1:]

      if len(initialVariables)>0:
        for i, var in enumerate(initialVariables):
          if "#" in var:
            continue
          v=json.loads(var)
          ivars.append(str(v[0]))
      if len(additionalVariables)>0:
        for i, var in enumerate(additionalVariables):
          if "#" in var:
            continue
          v=json.loads(var)
          avars.append(str(v[0]))
      if len(SpectatorVariables)>0:
        for i, var in enumerate(SpectatorVariables):
          if "#" in var:
            continue
          v=json.loads(var)
          svars.append(str(v[0]))
      #check the input variables
      print "initial Variables ", len(ivars), ivars
      print "additional Variables ", len(avars), avars
      print "spectator Variables", len(svars), svars
      allVars=avars+ivars+svars
      for var1 in allVars:
        counter=0
        for var2 in allVars:
          if var1==var2:
            counter+=1
          if counter>=2:
            print "ERROR: input variables defined multiple times"
            print var1
            exit(1)
      print "ivars ", ivars
      self.usedVariables=ivars
      self.unusedVariables=avars
      #self.TenBestMVAs[0]=(0.0,0.0,0.0,0.0,0.0,0,2,self.usedVariables, self.unusedVariables)
      #done with variables
      #read the rest of the stuff
      for line in configLines:
        if "#" in line:
          continue
        if "RunOn" in line:
          RunSystem=line.split("=",1)[1]
        if "NParticles" in line:
          self.nParticles=int(line.split("=",1)[1])
        if "NIterations" in line:
          self.NIterations=int(line.split("=",1)[1])
        if "wIneratia" in line:
          self.vw=float(line.split("=",1)[1])
        if "wMemory" in line:
          self.vp=float(line.split("=",1)[1])
        if "wSwarm" in line:
          self.vg=float(line.split("=",1)[1])
        if "FOM=" in line:
          self.FOM=line.split("=",1)[1]
        if "ADThreshold" in line:
          self.ADThreshold=float(line.split("=",1)[1])
        if "KSThreshold" in line:
          self.KSThreshold=float(line.split("=",1)[1])
        if "FindBestVariables" in line:
          self.FindBestVariables=int(line.split("=",1)[1])
        if "MaxVariablesInCombination" in line:
          self.MaxVariablesInCombination=int(line.split("=",1)[1])
        if "ImprovementThreshold" in line:
          self.ImprovementThreshold=float(line.split("=",1)[1])
        if "FactoryString" in line:
          self.FactoryString=line.split("=",1)[1]
        if "PreparationString" in line:
          self.PreparationString=line.split("=",1)[1]
        if "nTrain_B" in line:
          buff=line.split("=",1)[1]
          if buff!="" and buff!=" ":
            self.PreparationString="nTrain_B="+buff+":"+self.PreparationString
        if "nTrain_C" in line:
          buff=line.split("=",1)[1]
          if buff!="" and buff!=" ":
            self.PreparationString="nTrain_C="+buff+":"+self.PreparationString
        if "nTrain_O" in line:
          buff=line.split("=",1)[1]
          if buff!="" and buff!=" ":
            self.PreparationString="nTrain_O="+buff+":"+self.PreparationString
        if "nTest_C" in line:
          buff=line.split("=",1)[1]
          if buff!="" and buff!=" ":
            self.PreparationString="nTest_B="+buff+":"+self.PreparationString
        if "nTest_C" in line:
          buff=line.split("=",1)[1]
          if buff!="" and buff!=" ":
            self.PreparationString="nTest_C="+buff+":"+self.PreparationString
        if "nTest_O" in line:
          buff=line.split("=",1)[1]
          if buff!="" and buff!=" ":
            self.PreparationString="nTest_O="+buff+":"+self.PreparationString
        if "BWeightExpression" in line:
          self.BWeightExpression=line.split("=",1)[1]
        if "CWeightExpression" in line:
          self.CWeightExpression=line.split("=",1)[1]
        if "OWeightExpression" in line:
          self.OWeightExpression=line.split("=",1)[1]
          
        if "MethodType" in line:
          self.MethodType=line.split("=",1)[1]

        if "MethodParams" in line:
          self.MethodParams=line.split("=",1)[1]

        if 'UseFixedTrainTestSplitting' in line:

           arg = line.split('=', 1)[0].replace(' ','')

           if arg == 'UseFixedTrainTestSplitting':
              self.UseFixedTrainTestSplitting = line.split('=', 1)[1]

        if 'UseFixedTrainTestSplitting_Train' in line:

           arg = line.split('=', 1)[0].replace(' ','')

           if arg == 'UseFixedTrainTestSplitting_Train':
              self.UseFixedTrainTestSplitting_Train = line.split('=', 1)[1]

        if "SourceBTree" in line:
          self.BTreeName=line.split("=",1)[1]
        if "SourceCTree" in line:
          self.CTreeName=line.split("=",1)[1]
        if "SourceOTree" in line:
          self.OTreeName=line.split("=",1)[1]
        if "RepeatTrainingNTimes" in line:
          self.RepeatTrainingNTimes=int(line.split("=",1)[1])
        if "DrawNRandomAsStartingVars" in line:
          self.DrawNRandomAsStartingVars=int(line.split("=",1)[1])
        if "SaveTrainingsToTrees" in line:
          self.SaveTrainingsToTrees=line.split("=",1)[1]
        if "coord" in line:
          buff=line.split("=",1)[1]
          coord=json.loads(buff)
          print "setting up coordinate ",coord
          if coord not in self.Coordinates:
            self.Coordinates.append(coord)
            self.BestCoordinatesGlobal.append([coord[0],0.0])

      self.QueHelper = QueHelper(RunSystem)

    def InitParticles(self):

      subdir_idx_format = '0'+str(1+int(math.log10(self.nParticles)))+'d'

      for part in range(self.nParticles):

        part_dir = self.OutputDir+'/Particle'+format(part, subdir_idx_format)+'/'

        subprocess.call(['cp', '-r', self.DataSubdir, part_dir])
#        subprocess.call(['mv', part_dir+"/PSO.sh", part_dir+"/PSO"+format(part, subdir_idx_format)+".sh" ])

        #get starting values for the particle
        #uniformly distributed in coord space
        initialCoords=[]
        for coord in self.Coordinates:
          if coord[4]=="int":
            initValue=int(self.rand.Uniform(coord[1],coord[2]))
            initVelocity=int(self.rand.Uniform(-(coord[2]-coord[1]),(coord[2]-coord[1])))
            #now bound the velocity by the maximum
            initVelocity=int(cmp(initVelocity,0)*min(abs(initVelocity),coord[3]))
          elif coord[4]=="float":
            initValue=float(self.rand.Uniform(coord[1],coord[2]))
            initVelocity=float(self.rand.Uniform(-(coord[2]-coord[1]),(coord[2]-coord[1])))
            #now bound the velocity by the maximum
            initVelocity=float(cmp(initVelocity,0.0)*min(abs(initVelocity),coord[3]))
          else:
            print "ERROR: coordinates have to be int or float"
            exit(1)
          initialCoords.append([coord[0],initValue,initVelocity]) 
        print "Particle", part, " has inital coords ", initialCoords
        particle = Particle(
          part_dir,
          part,
          self.Verbose,
          self.usedVariables,
          self.unusedVariables,
          self.vw, self.vp, self.vg,
          self.Coordinates,
          initialCoords,
          self.FOM,
          self.ADThreshold,
          self.KSThreshold,
          self.FactoryString,
          self.PreparationString,
          self.BWeightExpression,
          self.CWeightExpression,
          self.OWeightExpression,
          self.BTreeName,
          self.CTreeName,
          self.OTreeName,
          self.MethodType,
          self.MethodParams,
          self.QueHelper,
          self.FindBestVariables,
          self.MaxVariablesInCombination,
          self.ImprovementThreshold,
          self.RepeatTrainingNTimes,
          self.DrawNRandomAsStartingVars,
          self.SaveTrainingsToTrees,
          self.UseFixedTrainTestSplitting,
          self.UseFixedTrainTestSplitting_Train
        )
#        particle.SetTestPoint(initTree,initShrinkage,initBagging,initCuts,2,1,0)
        self.Particles.append(particle)

#      print str(self.nParticles)+" Particles set up"

    def RunPSO(self):
      print "\n-------------------------------------------------------------"
      print "Starting Optimization"
      print "doing ", self.NIterations, " Iterations\n" 
      startTime=0
      finishTime=0
      totalTime=0.0
      nIterations=self.NIterations
      for it in range(nIterations):

        startTime = time.time()

        for particle in self.Particles:
            particle.StartEvaluation()

        running = True

        print "\nIteration ", it

        if it!=0:
           estTime=totalTime/it * (nIterations-it)/60.0/60.0
           print "optimization finished in ca ", estTime, " hours"

        check_dt_sec = int(60)

        print 'Number of particles finished: [checked every '+str(check_dt_sec)+'sec]'

        while running:

          nFinished = 0

          time.sleep(check_dt_sec)

          htc_jobIDs = None

          while htc_jobIDs == None:

             htc_jobIDs = HTCondor_jobIDs(os.environ['USER'], permissive=True)

             if htc_jobIDs == None:

                WARNING('PSOManager.py -- call to "condor_q" failed, will wait 60sec and try again')

                time.sleep(60)

          for i_part in self.Particles:

              i_part.ManageJob(htc_jobIDs)

              if not i_part.isRunning: nFinished += 1

          #print nFinished
          sys.stdout.write(str(nFinished)+" ")
          sys.stdout.flush()
          #print(str(nFinished)+" "),

          if nFinished == self.nParticles:
             running = False

        print " "
        for particle in self.Particles:

          fom, ad, ks, methodString, currentcoords, usedVars, unusedVars = particle.GetResult()

          if self.Verbose:
             print "particle returned: "
             print fom, ad, ks, methodString, currentcoords, usedVars, unusedVars

          if fom >= self.TenBestMVAs[0][0]:

             if fom==0.0:
                self.TenBestMVAs[0] = [fom, ad, ks, methodString, usedVars, unusedVars,currentcoords]

             else:
                self.TenBestMVAs[9] = [fom, ad, ks, methodString, usedVars, unusedVars,currentcoords]
                self.TenBestMVAs = sorted(self.TenBestMVAs, key=lambda s:s[0], reverse=True)

          elif fom >= self.TenBestMVAs[9][0]:
             self.TenBestMVAs[9] = [fom, ad, ks, methodString, usedVars, unusedVars,currentcoords]
             self.TenBestMVAs = sorted(self.TenBestMVAs, key=lambda s:s[0], reverse=True)

        self.BestFOMGlobal         = self.TenBestMVAs[0][0]
        self.BestADGlobal          = self.TenBestMVAs[0][1]
        self.BestKSGlobal          = self.TenBestMVAs[0][2]
        self.BestCoordinatesGlobal = self.TenBestMVAs[0][6]

        for particle in self.Particles:
            particle.UpdateParticle(self.BestCoordinatesGlobal, it, self.BestFOMGlobal, self.BestADGlobal, self.BestKSGlobal)

        print "\n------------------------------------------------------------------------"

        print "Best Result after Iteration "+str(it)

        print self.TenBestMVAs[0][:6]

        #print "testing blabla "
        print "OutputDir = ", self.OutputDir

        self.SaveStatus(self.OutputDir+"/PSOResult.txt", self.OutputDir+"/FinalMVAConfig_PSO.txt", self.OutputDir+'.conf')

        finishTime = time.time()

        totalTime += (finishTime-startTime)

    def PrintResult(self):
        print 'Ten Best MVAs'
        for i in range(10):
            print self.TenBestMVAs[i][:6]

    def SaveStatus(self, SaveFile, BestBDTFile, FinalMVAConfFile=None):
        savefile = open(SaveFile, 'w')
        savefile.write('nParticles '     +str(self.nParticles)     +'\n')
        savefile.write('wIneratia '      +str(self.vw)             +'\n')
        savefile.write('wMemory '        +str(self.vp)             +'\n')
        savefile.write('wSwarm '         +str(self.vg)             +'\n')
        savefile.write('OutputDir '      +str(self.OutputDir)      +'\n')
        savefile.write('usedVariables '  +str(self.usedVariables)  +'\n')
        savefile.write('unusedVariables '+str(self.unusedVariables)+'\n')
        savefile.write('ADThreshold '    +str(self.ADThreshold)    +'\n')
        savefile.write('KSThreshold '    +str(self.KSThreshold)    +'\n')
        savefile.write('TenBestMVAs\n')
        for i in range(10): savefile.write(str(self.TenBestMVAs[i][:6])+'\n')
        savefile.write('\n')
        savefile.close()

#        for part in self.Particles: part.SaveParticleStatus()

        bestBDTFile = open(BestBDTFile, 'w')
        bestBDTFile.write('FOM '         +str(self.TenBestMVAs[0][0])         +'\n')
        bestBDTFile.write('AD '          +str(self.TenBestMVAs[0][1])         +'\n')
        bestBDTFile.write('KS '          +str(self.TenBestMVAs[0][2])         +'\n')
        bestBDTFile.write('Method '      +str(self.TenBestMVAs[0][3])         +'\n')
        bestBDTFile.write('Variables '   +str(self.TenBestMVAs[0][4])         +'\n')
        bestBDTFile.write('Factory '     +str(self.FactoryString)             +'\n')
        bestBDTFile.write('Preparation ' +str(self.PreparationString)         +'\n')
        bestBDTFile.write('BWeight '+str(self.BWeightExpression)    +'\n')
        bestBDTFile.write('CWeight '+str(self.CWeightExpression)    +'\n')
        bestBDTFile.write('OWeight '+str(self.OWeightExpression)    +'\n')
        #bestBDTFile.write('Preparation ' +str(self.BackgroundWeightExpression)+'\n')
        bestBDTFile.write(                str(self.TenBestMVAs[0][:6])        +'\n')
        bestBDTFile.write('\n')
        bestBDTFile.close()

        print " OutputDir = " , self.OutputDir
        print " usedVariables = ", self.usedVariables

        if FinalMVAConfFile != None:

           # disable Silent mode in FactoryString of .conf file
           conf_FactoryString = str(self.FactoryString)

           if conf_FactoryString.startswith('Silent:'):
              conf_FactoryString = '!Silent:'+conf_FactoryString[len(':Silent'):]

           if conf_FactoryString.endswith(':Silent'):
              conf_FactoryString = conf_FactoryString[:-len(':Silent')]+':!Silent'

           if ':Silent:' in conf_FactoryString:
              conf_FactoryString = conf_FactoryString.replace(':Silent:', ':!Silent:')

           finalMVAConfFile = open(FinalMVAConfFile, 'w')
           finalMVAConfFile.write('[configuration]'+'\n')
           finalMVAConfFile.write('\n')
           finalMVAConfFile.write('factory    = '+str(conf_FactoryString)             +'\n')
           finalMVAConfFile.write('dataloader = '+str(self.PreparationString)         +'\n')
           finalMVAConfFile.write('\n')
           finalMVAConfFile.write('method = kBDT, BDTG, '+str(self.TenBestMVAs[0][3]) +'\n')
           finalMVAConfFile.write('\n')
           finalMVAConfFile.write('[variables]'+'\n')
           for i_var in self.TenBestMVAs[0][4]: finalMVAConfFile.write(str(i_var)+'\n')
           finalMVAConfFile.write('\n')
           finalMVAConfFile.write('[spectators]'+'\n')
           for i_spc in self.TenBestMVAs[0][5]: finalMVAConfFile.write(str(i_spc)+'\n')
           finalMVAConfFile.close()

    def GetVariableNumbers(self):
      outfile=open(self.OutputDir+"/VariableNumbers.txt","w")
      for var in self.usedVariables:
        VarCount=0
        for i in range(self.nParticles):
          file=open(self.OutputDir+"/Particle"+str(i)+"/ParticleRoute.txt","r")
          lines=list(file)
          for line in lines:
            #print line
            if "0.0 " not in line:
              #print line
              if var in line:
                #print line
                VarCount+=1
          file.close()
        print var, VarCount
        outfile.write(var+" "+str(VarCount)+"\n")
      for var in self.unusedVariables:
        VarCount=0
        for i in range(self.nParticles):
          file=open(self.OutputDir+"/Particle"+str(i)+"/ParticleRoute.txt","r")
          lines=list(file)
          for line in lines:
            if "0.0 " not in line:
              if var in line:
                VarCount+=1
          file.close()
        print var, VarCount
        outfile.write(var+" "+str(VarCount)+"\n")
      outfile.close()

#    def testFunction(self):
#      for part in self.Particles:
#        part.WriteConfig()
#        part.UpdateParticle([[u'NTrees', 1000], [u'Shrinkage', 0.02], [u'BaggedSampleFraction', 0.2], [u'nCuts', 50]])
#        part.StartEvaluation()
#        part.CheckJobStatus()
#        print self.BestCoordinatesGlobal
#        bg=[[u'NTrees', 0.0], [u'Shrinkage', 0.0], [u'BaggedSampleFraction', 0.0], [u'nCuts', 0.0]]
#
#      print "did the tests"

    def CompileAndSetupClientExecutable(self):

        if not os.path.isdir(self.DataSubdir):
           raise RuntimeError('target directory for Particle executable not found: '+self.DataSubdir)

        Particle_C = 'PSO/Particle.C'
        if not os.path.isfile(Particle_C):
           raise RuntimeError('source code not found: '+Particle_C)

        print 'Compiling '+Particle_C
        print('g++ -o '+self.DataSubdir+'/Particle '+Particle_C+' `root-config --cflags --glibs` -lTMVA')
        subprocess.call(['g++ -o '+self.DataSubdir+'/Particle '+Particle_C+' `root-config --cflags --glibs` -lTMVA'], shell=True)
