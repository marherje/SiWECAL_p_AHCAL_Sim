#!/usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess, math, ROOT

from array import array

from common import *

class Particle:

    def __init__(self, Path,particleNumber,Verbose, usedVariables, unusedVariables, vw, vp, vg, coordinates, initialcoordinates, FOM, ADThreshold, KSThreshold, FactoryString, PreparationString, BWeightExpression, CWeightExpression, OWeightExpression, BTreeName, CTreeName, OTreeName, MethodType, MethodParams, QueHelper, FindBestVariables, MaxVariablesInCombination, ImprovementThreshold, RepeatTrainingNTimes, DrawNRandomAsStartingVars,SaveTrainingsToTrees,UseFixedTrainTestSplitting,UseFixedTrainTestSplitting_Train):

      self.particleNumber=particleNumber
      self.Iteration=0

      self.Path=Path

      self.fpath_confHTC     = None
      self.fpath_confHTC_err = None
      self.fpath_confHTC_out = None
      self.fpath_confHTC_log = None

      self.ExeFile = None
      self.RunFile = self.Path+'/run.sh'

      self.Verbose = Verbose
      self.JobID = None
      self.isRunning = False
      self.rand=ROOT.TRandom3(int(self.particleNumber)) ## Use seed=0 to get non reproduceability

      self.initialVariables=usedVariables
      self.additionalVariables=unusedVariables
      self.LastUsedVariables=[]
      self.LastUnusedVariables=[]
      self.vw=vw
      self.vp=vp
      self.vg=vg
      self.Coordinates=coordinates
      self.BestCoordinates=[]
      self.BestCoordinatesGlobal=[]

      for coord in self.Coordinates:
          self.BestCoordinates.append([coord[0],0.0])
          self.BestCoordinatesGlobal.append([coord[0],0.0])

      self.currentCoordinates=initialcoordinates
      self.FOM=FOM
      self.ADThreshold=ADThreshold
      self.KSThreshold=KSThreshold
      self.FactoryString=FactoryString
      self.PreparationString=PreparationString
      self.BWeightExpression=BWeightExpression
      self.CWeightExpression=CWeightExpression
      self.OWeightExpression=OWeightExpression
      self.BTreeName=BTreeName
      self.CTreeName=CTreeName
      self.OTreeName=OTreeName
      self.MethodType=MethodType
      self.MethodParams=MethodParams
      self.QueHelper=QueHelper
      self.FindBestVariables=FindBestVariables
      self.MaxVariablesInCombination=MaxVariablesInCombination
      self.ImprovementThreshold=ImprovementThreshold
      self.RepeatTrainingNTimes=RepeatTrainingNTimes
      self.DrawNRandomAsStartingVars=DrawNRandomAsStartingVars

      if self.Verbose:
         print "-------------------------------------------"
         print "setting up particle", self.particleNumber

      if self.DrawNRandomAsStartingVars>0:

         allVars=self.initialVariables+self.additionalVariables

         il=len(allVars)

         self.initialVariables=[]

         for i in range(self.DrawNRandomAsStartingVars):
             idx=int(self.rand.Uniform(0,len(allVars)))
             var=allVars.pop(idx)
             self.initialVariables.append(var)

         self.additionalVariables=allVars

         if self.Verbose:
            print "Drawing ", self.DrawNRandomAsStartingVars, " starting Variables"
            print "usedVars\n", self.initialVariables
            print "unusedVariables\n", self.additionalVariables
            print il, len(self.initialVariables), len(self.additionalVariables)

      self.SaveTrainingsToTrees = SaveTrainingsToTrees

      self.UseFixedTrainTestSplitting       = UseFixedTrainTestSplitting
      self.UseFixedTrainTestSplitting_Train = UseFixedTrainTestSplitting_Train

      self.AllVariablesAtStart=self.additionalVariables+self.initialVariables
      self.AllVariablesAfterIteration=self.AllVariablesAtStart

      self.BestFOM=0.0
      self.BestAD=0.0
      self.BestKS=0.0
      self.BestFOMGlobal=0.0
      self.BestADGlobal=0.0
      self.BestKSGlobal=0.0

      # write execution script and run script for job submission
      self.ExeFile = self.Path+'/PSO'+str(self.particleNumber)+'.sh'

      exe_file = open(self.ExeFile, 'w')

      for line in self.QueHelper.GetExecLines():
          exe_file.write(line)

      exe_file.write('cd '+self.Path+'\n')
      exe_file.write('./Particle'   +'\n')

      exe_file.close()

      subprocess.call(['chmod', 'u+x', self.ExeFile])

      if self.QueHelper.GetConfigLines() == []:

         run_file = open(self.RunFile, 'w')

         for line in self.QueHelper.GetRunLines():
             line = line.replace("INSERTPATHHERE"      , self.Path)
             line = line.replace("INSERTEXECSCRIPTHERE", self.Path+"/PSO"+str(self.particleNumber)+".sh")
             run_file.write(line)

         run_file.close()

         subprocess.call(['chmod', 'u+x', self.RunFile])

      else:

         # HTCondor configuration script
         self.fpath_confHTC = os.path.abspath(self.Path+'/conf.htc')

         config_file = open(self.fpath_confHTC, 'w')
         config_lines = self.QueHelper.GetConfigLines()

         for line in config_lines:

             line = line.replace('__PATH__'      , self.Path)
             line = line.replace('__EXEC_FILE__' , self.Path+'/PSO'+str(self.particleNumber)+'.sh')

             line = line.replace('__NAME__'      , 'PSO'+str(self.particleNumber))
             line = line.replace('__BATCH_NAME__', 'PSO'+str(self.particleNumber))

             config_file.write('\n'+line+'\n')

         config_file.close()

         # execution script
         run_file = open(self.RunFile, 'w')

         run_file.write('#!/bin/sh')
         run_file.write('\n\n')
         run_file.write('cd '+os.path.dirname(self.fpath_confHTC))
         run_file.write('\n')
         run_file.write('condor_submit '+os.path.basename(self.fpath_confHTC))
         run_file.write('\n')

         run_file.close()

         subprocess.call(['chmod', 'u+x', self.RunFile])

      self.WriteConfig()

    def WriteConfig(self):
      # write ConfigFile
      configfile = open(self.Path+"/ParticleConfig.txt","w")
      configfile.write("particleNumber "+str(self.particleNumber)+"\n")
      configfile.write("Iteration "+str(self.Iteration)+"\n")
      configfile.write("FOM "+str(self.FOM)+"\n")
      configfile.write("SaveTrainingsToTrees "+str(self.SaveTrainingsToTrees)+"\n")
      configfile.write("ADThreshold "+str(self.ADThreshold)+"\n")
      configfile.write("KSThreshold "+str(self.KSThreshold)+"\n")
      configfile.write("FactoryString "+str(self.FactoryString)+"\n")
      configfile.write("PreparationString "+str(self.PreparationString)+"\n")
      configfile.write("BWeightExpression "+str(self.BWeightExpression)+"\n")
      configfile.write("CWeightExpression "+str(self.CWeightExpression)+"\n")
      configfile.write("OWeightExpression "+str(self.OWeightExpression)+"\n")
      configfile.write("BTreeName "+str(self.BTreeName)+"\n")
      configfile.write("CTreeName "+str(self.CTreeName)+"\n")
      configfile.write("OTreeName "+str(self.OTreeName)+"\n")

      configfile.write("UseFixedTrainTestSplitting "+str(self.UseFixedTrainTestSplitting)+"\n")

      str_UseFixedTrainTestSplitting_Train = str(self.UseFixedTrainTestSplitting_Train)

      while str_UseFixedTrainTestSplitting_Train.startswith('"') or str_UseFixedTrainTestSplitting_Train.startswith("'"):
            str_UseFixedTrainTestSplitting_Train = str_UseFixedTrainTestSplitting_Train[1:]

      while str_UseFixedTrainTestSplitting_Train.endswith('"') or str_UseFixedTrainTestSplitting_Train.endswith("'"):
            str_UseFixedTrainTestSplitting_Train = str_UseFixedTrainTestSplitting_Train[:-1]

      configfile.write("UseFixedTrainTestSplitting_Train \""+str_UseFixedTrainTestSplitting_Train+"\"\n")

      configfile.write("FindBestVariables "+str(self.FindBestVariables)+"\n")
      configfile.write("MaxVariablesInCombination "+str(self.MaxVariablesInCombination)+"\n")
      configfile.write("ImprovementThreshold "+str(self.ImprovementThreshold)+"\n")
      configfile.write("RepeatTrainingNTimes "+str(self.RepeatTrainingNTimes)+"\n")
      configfile.write("MethodType "+str(self.MethodType)+"\n")
      methodString=self.MethodParams
      # print methodString
      for coord in self.currentCoordinates:
        configfile.write("coord "+coord[0]+" "+str(coord[1])+"\n")
        if coord[0] in methodString:
          splitMethods=methodString.split(":")
          for i,par in enumerate(splitMethods):
            if par.split("=",1)[0]==coord[0]:
              val=par.split("=",1)[1]
              splitMethods[i]=par.replace(val,str(coord[1]))
              break
          newMethods=":"
          methodString=newMethods.join(splitMethods)
      self.MethodParams=methodString
      configfile.write("MethodParameters "+str(self.MethodParams)+"\n")

      if len(self.initialVariables) == 0:
         KILL('particle.py -- WriteConfig: empty list of initial variables')

      configfile.write("--InitialVariables--"+"\n")
      for var in self.initialVariables: configfile.write(var+"\n")
      configfile.write("--EndInitVars--\n")

      configfile.write("--AdditionalVariables--"+"\n")
      for var in self.additionalVariables: configfile.write(var+"\n")
      configfile.write("--EndAddVars--\n")

      configfile.close()

    def StartEvaluation(self):

        self.JobID = self.QueHelper.StartJob(self.Path+'/run.sh')
#        print self.JobID

        fpath_dc = {}

        with open(self.fpath_confHTC, 'r') as cfg_file:

             cfg_lines = cfg_file.readlines()

             for i_line in cfg_lines:
                 i_line = i_line.replace('\n', '')
                 i_line = i_line.replace(' ' , '')

                 if   i_line.startswith('error=') : fpath_dc['err'] = i_line[len('error='):]
                 elif i_line.startswith('output='): fpath_dc['out'] = i_line[len('output='):]
                 elif i_line.startswith('log=')   : fpath_dc['log'] = i_line[len('log='):]

        for i_type in ['err', 'out', 'log']:

            if i_type not in fpath_dc:
               KILL('particle.py -- StartEvaluation: path to job \"'+i_type+'\" file not found in HTCondor configuration file: '+self.fpath_confHTC)

        self.fpath_confHTC_err = fpath_dc['err'].replace('$(Cluster).$(Process)', str(self.JobID))
        self.fpath_confHTC_out = fpath_dc['out'].replace('$(Cluster).$(Process)', str(self.JobID))
        self.fpath_confHTC_log = fpath_dc['log'].replace('$(Cluster).$(Process)', str(self.JobID))

        self.isRunning = True

        return

    def ManageJob(self, jobID_dict):

        # delete_test = ['Error in <TH1F::Smooth>: Smooth only supported for histograms with >= 3 bins. Nbins = 2','Error in <TH1F::Smooth>: Smooth only supported for histograms with >= 3 bins. Nbins = 1','Error in <TDecompLU::InvertLU>: matrix is singular, 1 diag elements < tolerance of 2.2204e-16']
        # Patch to remove the smoothing bug START
        if os.path.isfile(self.fpath_confHTC_err) and (os.path.getsize(self.fpath_confHTC_err) != 0):
            delete_test = ['Error in <TH1F::Smooth>: Smooth only supported for histograms with >= 3 bins. Nbins = 2','Error in <TH1F::Smooth>: Smooth only supported for histograms with >= 3 bins. Nbins = 1']
            f = open(str(self.fpath_confHTC_err),'r')
            lst = []
            for line in f:
                for string in delete_test:
                    if string in line:
                        line = line.replace(string,'')
                        if line != '\n':
                            lst.append(line)
            f.close()
            f = open(str(self.fpath_confHTC_err),'w')
            for line in lst:
                f.write(line)
            f.close()
        # Patch to remove the smoothing bug END
        # Change the condition of !=0 to > 100 so the error is at least one sentence and

        if os.path.isfile(self.fpath_confHTC_err) and (os.path.getsize(self.fpath_confHTC_err) != 0):

           log_msg = 'job '+str(self.JobID)+' returned ERROR (stderr file not empty),'

           self.QueHelper.KillJob(str(self.JobID))

           # resubmit
           self.StartEvaluation()

           log_msg += ' has been removed and resubmitted as job '+str(self.JobID)

           WARNING('particle.py -- ManageJob: '+log_msg)

        else:

           job_exists = bool(str(self.JobID) in jobID_dict)

           if not job_exists:

              self.isRunning = False

           else:

              job_status = jobID_dict[str(self.JobID)]['STATUS']

              if bool(job_status == 'C'):

                 self.isRunning = False

              elif bool(job_status in ['H','X']):

                 log_msg = 'job '+str(self.JobID)+' found in status "'+job_status+'"'

                 self.QueHelper.KillJob(str(self.JobID))

                 # resubmit
                 self.StartEvaluation()

                 log_msg += ' has been removed and resubmitted as job '+str(self.JobID)

                 WARNING('particle.py -- ManageJob: '+log_msg)

              else:

                 self.isRunning = True

        return

    def UpdateParticle(self, BestCoordsGlobal,Iteration,bestFOMGlobal, bestADGlobal, bestKSGlobal):
      self.Iteration=Iteration
      self.BestCoordinatesGlobal=BestCoordsGlobal
      self.BestFOMGlobal=bestFOMGlobal
      self.BestADGlobal=bestADGlobal
      self.BestKSGlobal=bestKSGlobal

      if self.Verbose:
         print "\nUpdating particle ",self.particleNumber
         print "current ", self.currentCoordinates
         print "best global ", self.BestCoordinatesGlobal
         print "best particle ", self.BestCoordinates

      newCoords=[]
      for coord in self.Coordinates:
        curVel=0
        newVel=0
        cc="bla"
        newcoord=""
        bcg=""
        bcp=""
        for c in self.BestCoordinatesGlobal:
          if coord[0]==c[0]:
            bcg=c[1]
            break
        for c in self.BestCoordinates:
          if coord[0]==c[0]:
            bcp=c[1]
            break
        for c in self.currentCoordinates:
          if coord[0]==c[0]:
            cc=c[1]
            curVel=c[2]
            break

        rp=self.rand.Rndm()
        rg=self.rand.Rndm()

        newVel = self.vw*curVel + self.vp*rp*(bcp-cc)+self.vg*rg*(bcg-cc)

        newVel = cmp(newVel,0) * min(abs(newVel),coord[3])

        if   coord[4] == 'int'  : newVel = int  (newVel)
        elif coord[4] == 'float': newVel = float(newVel)

        newcoord = cc + newVel
        newcoord = abs(cmp(newcoord,0)*min(abs(newcoord),coord[2]))
        newcoord = abs(cmp(newcoord,0)*max(abs(newcoord),coord[1]))

        if self.Verbose:
           print "\nOld Coordinate ", coord[0],cc,curVel
           print "best global ", bcg
           print "best for this particle ", bcp

           print "New Coordinate ", newcoord, newVel

        newCoords.append([coord[0],newcoord,newVel])

      self.currentCoordinates=newCoords

      self.WriteConfig()

    def GetResult(self):
      resultfile = open(self.Path+"/ParticleResult.txt","r")
      lines = list(resultfile)

      fom = 0.0
      ad  = 0.0
      ks = 0.0

      self.LastUsedVariables=self.initialVariables
      self.LastUnusedVariables=self.additionalVariables
      self.initialVariables=[]
      self.additionalVariables=[]

      for line in lines:
        if "BestFOM" in line:
          fom=float(line.split(" ",1)[1])
        if "ADScore" in line:
          ad=float(line.split(" ",1)[1])
        if "KSScore" in line:
          ks=float(line.split(" ",1)[1])
        if "MethodString" in line:
          self.MethodParams=line.split(" ",1)[1].strip()
        if "UsedVar" in line:
          self.initialVariables.append(line.split(" ",1)[1].strip())
        if "UnusedVars" in line:
          self.additionalVariables.append(line.split(" ",1)[1].strip())

      resultfile.close()

      #self.AllVariablesAfterIteration=self.initialVariables+self.additionalVariables
      self.additionalVariables=[]
      for vv in self.AllVariablesAtStart:
        if vv not in self.initialVariables:
          self.additionalVariables.append(vv)

      if ad<self.ADThreshold or ks<self.KSThreshold:
         fom = 0.0
         ad  = 0.0
         ks = 0.0

      else:
        if fom>=self.BestFOM:

           self.BestFOM=fom
           self.BestAD=ad
           self.BestKS=ks
           self.BestCoordinates=[]

           for coord in self.currentCoordinates:
               self.BestCoordinates.append([coord[0],coord[1]])

      if self.Verbose:
         print "particle ", self.particleNumber
         print self.BestCoordinates
         print self.BestFOM

      RouteFile=open(self.Path+"/ParticleRoute.txt","a")
      Route=str(fom).replace("\n","")+" "+str(ad).replace("\n","")+" "+str(ad).replace("\n","")+" "
      for ccc in self.currentCoordinates:
        Route+=str(ccc[1])+" "
      Route+=str(self.additionalVariables)+"\n"
      RouteFile.write(Route)
      RouteFile.write("--Next--\n")
      RouteFile.close()

      return fom, ad, ks, self.MethodParams, self.currentCoordinates, self.initialVariables, self.additionalVariables

#    def SaveParticleStatus(self):
#        PartSavefile=open(self.Path+"/ParticleStatus.txt","w")
##        PartSavefile.write(str(self.Path)+"\n")
#        PartSavefile.write(str(self.KSThreshold)+"\n")
#        PartSavefile.write(str(self.usedVariables)+"\n")
#        PartSavefile.write(str(self.unusedVariables)+"\n")
#        PartSavefile.write(str(self.particleNumber)+"\n")
#        PartSavefile.write(str(self.BestROC)+"\n")
#        PartSavefile.write(str(self.BestKS)+"\n")
#        PartSavefile.write(str(self.BestNTrees)+"\n")
#        PartSavefile.write(str(self.BestShrinkage)+"\n")
#        PartSavefile.write(str(self.BestBagging)+"\n")
#        PartSavefile.write(str(self.BestCuts)+"\n")
#        PartSavefile.write(str(self.BestDepth)+"\n")
#        PartSavefile.write(str(self.VelTree)+"\n")
#        PartSavefile.write(str(self.VelShrinkage)+"\n")
#        PartSavefile.write(str(self.VelBagging)+"\n")
#        PartSavefile.write(str(self.VelCuts)+"\n")
#        PartSavefile.write(str(self.VelDepth)+"\n")
#        PartSavefile.write(str(self.NTrees)+"\n")
#        PartSavefile.write(str(self.Shrinkage)+"\n")
#        PartSavefile.write(str(self.Bagging)+"\n")
#        PartSavefile.write(str(self.Cuts)+"\n")
#        PartSavefile.write(str(self.Depth)+"\n")
#        PartSavefile.close()
