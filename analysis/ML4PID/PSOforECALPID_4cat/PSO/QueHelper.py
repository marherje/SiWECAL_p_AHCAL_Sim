#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys, subprocess, time

from common import *

class QueHelper:
  def __init__(self,RunSystem):
    self.RunSystem=RunSystem

    self.ExecLines = []

    self.ConfigLines = []

    self.RunLines = []

    # change if you want to use a different CMSSW Version
    self.CMSSW_BASE = os.environ['PATH']
    self.SCRAM_ARCH = os.environ['PATH']
    print("\n----------------------------------------\n")
    print "using CMSSW_BASE="+self.CMSSW_BASE
    print "using SCRAM_ARCH="+self.SCRAM_ARCH
    print "might not be the right choice, depending on your linux, target linux and CMSSW Version"
    print "changeable in QueHelper.py"

    if RunSystem=="EKPSL5":
      thisPortal=os.environ["HOSTNAME"]
      if thisPortal=="ekpcms5":
        thisque=os.environ["SGE_CLUSTER_NAME"]
        if thisque!="p_ogs1111":
          print "you need to setup the sl5 que first"
          print "source /opt/sge62/ekpcluster/common/settings.sh"
          exit(1)
      elif thisPortal=="ekpcms6":
        print "using sl5 que on ekpcms6 might lead to problems"
        thisque=os.environ["SGE_CLUSTER_NAME"]
        if thisque!="p_ogs1111":
          print "you need to setup the sl5 que first"
          print "source /opt/sge62/ekpcluster/common/settings.sh"
          exit(1)
      else:
        print "dont know this portal"
        exit(1)
      self.ExecLines=[
          "export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n",
          "export SCRAM_ARCH="+self.SCRAM_ARCH+"\n",
          "source $VO_CMS_SW_DIR/cmsset_default.sh\n",
          "cd "+self.CMSSW_BASE+"/src\n",
          "eval `scram runtime -sh`\n"
          ]
      self.RunLines=[
          "qsub -cwd -S /bin/bash -o INSERTPATHHERE/logs/\$JOB_NAME.o\$JOB_ID -e INSERTPATHHERE/logs/\$JOB_NAME.e\$JOB_ID -q 'medium' INSERTEXECSCRIPTHERE\n"
          ] 
    elif RunSystem=="EKPSL6":
      thisPortal=os.environ["HOSTNAME"]
      #print thisPortal
      if thisPortal!="ekpcms6":
        print "WARNING"
        print "you try to start jobs on the ekp SL6 que from ekpcms5"
        print "do it from ekpcms or manually change QueHelper.py"
      thisque=os.environ["SGE_CLUSTER_NAME"]
      if thisque!="p_ogs1111_sl6":
          print "you need to setup the sl6 que first"
          print "source /opt/ogs_sl6/ekpclusterSL6/common/settings.sh"
          exit(1)
      self.ExecLines=[
          "export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n",
          "export SCRAM_ARCH="+self.SCRAM_ARCH+"\n",
          "source $VO_CMS_SW_DIR/cmsset_default.sh\n",
          "cd "+self.CMSSW_BASE+"/src\n",
          "eval `scram runtime -sh`\n"
          ]
      self.RunLines=[
          "qsub -cwd -S /bin/bash -o INSERTPATHHERE/logs/\$JOB_NAME.o\$JOB_ID -e INSERTPATHHERE/logs/\$JOB_NAME.e\$JOB_ID -q 'medium' INSERTEXECSCRIPTHERE\n"
          ]
    elif RunSystem == "NAFSL7":

      self.ExecLines = [

        "#!/bin/bash\n",
        "source /etc/profile.d/modules.sh\n",
        "module use -a /afs/desy.de/group/cms/modulefiles/\n",
        "#module load cmssw/"+self.SCRAM_ARCH+"\n",
        "export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n",
        "export SCRAM_ARCH="+self.SCRAM_ARCH+"\n",
        "source $VO_CMS_SW_DIR/cmsset_default.sh\n",
        "cd "+self.CMSSW_BASE+"/src\n",
        "eval `scram runtime -sh`\n",
      ]

      # HTCondor getenv=True does not export LD_LIBRARY_PATH
      # --> added by hand in the script itself
      if 'LD_LIBRARY_PATH' in os.environ:
         self.ExecLines += ['\n']
         self.ExecLines += ['export LD_LIBRARY_PATH='+os.environ['LD_LIBRARY_PATH']]
         self.ExecLines += ['\n\n']

      self.ConfigLines = [

        'batch_name = __BATCH_NAME__',

        'executable = __EXEC_FILE__',

        'output = __PATH__/logs/__NAME__'+'.out.'+'$(Cluster).$(Process)',
        'error  = __PATH__/logs/__NAME__'+'.err.'+'$(Cluster).$(Process)',
        'log    = __PATH__/logs/__NAME__'+'.log.'+'$(Cluster).$(Process)',

        '#arguments = ',

        'transfer_executable = True',

        'universe = vanilla',

        'getenv = True',

        'should_transfer_files   = IF_NEEDED',
        'when_to_transfer_output = ON_EXIT',

        '#requirements = (OpSysAndVer == "SL6")',
        'requirements = (OpSysAndVer == "CentOS7")',

        ' RequestMemory  =  2000',
        '+RequestRuntime = 10799',

        'queue',
      ]

    elif RunSystem == "NAFSL6":

      self.ExecLines = [

        "#!/bin/bash\n",
        "source /etc/profile.d/modules.sh\n",
        "module use -a /afs/desy.de/group/cms/modulefiles/\n",
        "#module load cmssw/"+self.SCRAM_ARCH+"\n",
        "export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n",
        "export SCRAM_ARCH="+self.SCRAM_ARCH+"\n",
        "source $VO_CMS_SW_DIR/cmsset_default.sh\n",
        "cd "+self.CMSSW_BASE+"/src\n",
        "eval `scram runtime -sh`\n",
      ]

      # HTCondor getenv=True does not export LD_LIBRARY_PATH
      # --> added by hand in the script itself
      if 'LD_LIBRARY_PATH' in os.environ:
         self.ExecLines += ['\n']
         self.ExecLines += ['export LD_LIBRARY_PATH='+os.environ['LD_LIBRARY_PATH']]
         self.ExecLines += ['\n\n']

      self.ConfigLines = [

        'batch_name = __BATCH_NAME__',

        'executable = __EXEC_FILE__',

        'output = __PATH__/logs/__NAME__'+'.out.'+'$(Cluster).$(Process)',
        'error  = __PATH__/logs/__NAME__'+'.err.'+'$(Cluster).$(Process)',
        'log    = __PATH__/logs/__NAME__'+'.log.'+'$(Cluster).$(Process)',

        '#arguments = ',

        'transfer_executable = True',

        'universe = vanilla',

        'getenv = True',

        'should_transfer_files   = IF_NEEDED',
        'when_to_transfer_output = ON_EXIT',

        'requirements = (OpSysAndVer == "SL6")',
        '#requirements = (OpSysAndVer == "CentOS7")',

        ' RequestMemory  =  2000',
        '+RequestRuntime = 10799',

        'queue',
      ]

    elif RunSystem == "NAFSL7":

        self.ExecLines = [

          "#!/bin/bash\n",
          "source /etc/profile.d/modules.sh\n",
          "module use -a /afs/desy.de/group/cms/modulefiles/\n",
          "#module load cmssw/"+self.SCRAM_ARCH+"\n",
          "export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n",
          "export SCRAM_ARCH="+self.SCRAM_ARCH+"\n",
          "source $VO_CMS_SW_DIR/cmsset_default.sh\n",
          "cd "+self.CMSSW_BASE+"/src\n",
          "eval `scram runtime -sh`\n",
        ]

        # HTCondor getenv=True does not export LD_LIBRARY_PATH
        # --> added by hand in the script itself
        if 'LD_LIBRARY_PATH' in os.environ:
           self.ExecLines += ['\n']
           self.ExecLines += ['export LD_LIBRARY_PATH='+os.environ['LD_LIBRARY_PATH']]
           self.ExecLines += ['\n\n']

        self.ConfigLines = [

          'batch_name = __BATCH_NAME__',

          'executable = __EXEC_FILE__',

          'output = __PATH__/logs/__NAME__'+'.out.'+'$(Cluster).$(Process)',
          'error  = __PATH__/logs/__NAME__'+'.err.'+'$(Cluster).$(Process)',
          'log    = __PATH__/logs/__NAME__'+'.log.'+'$(Cluster).$(Process)',

          '#arguments = ',

          'transfer_executable = True',

          'universe = vanilla',

          'getenv = True',

          'should_transfer_files   = IF_NEEDED',
          'when_to_transfer_output = ON_EXIT',

          '#requirements = (OpSysAndVer == "SL6")',
          'requirements = (OpSysAndVer == "CentOS7")',

          ' RequestMemory  =  2000',
          '+RequestRuntime = tomorrow',

          'queue',
        ]

    elif RunSystem=="NAFSL5":
      self.ExecLines=[
        "#!/bin/bash",
        ". /etc/profile.d/modules.sh\n",
        "module use -a /afs/desy.de/group/cms/modulefiles/\n",
        "module load cmssw/"+self.SCRAM_ARCH+"\n",
        "export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n",
        "export SCRAM_ARCH="+self.SCRAM_ARCH+"\n",
        "source $VO_CMS_SW_DIR/cmsset_default.sh\n",
        "cd "+self.CMSSW_BASE+"/src\n",
        "eval `scram runtime -sh`\n"
      ]

      qsub_opts = [
        '-S /bin/bash',
        '-V',
#        '-pe local 4-8',
        '-l h_vmem=4G',
        '-l h_fsize=1G',
        '-l h_rt=96:00:00',
      ]

      self.RunLines = [
        'qsub '+' '.join(qsub_opts)+' -l os=sld5 -o INSERTPATHHERE/logs/\$JOB_NAME.o\$JOB_ID -e INSERTPATHHERE/logs/\$JOB_NAME.e\$JOB_ID INSERTEXECSCRIPTHERE\n'
      ] 

    elif RunSystem == "LXPLUS":

      self.ExecLines = [

        "#!/bin/bash\n",
        "source /etc/profile.d/modules.sh\n",
        "module use -a /afs/desy.de/group/cms/modulefiles/\n",
        "#module load cmssw/"+self.SCRAM_ARCH+"\n",
        "export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n",
        "export SCRAM_ARCH="+self.SCRAM_ARCH+"\n",
        "source $VO_CMS_SW_DIR/cmsset_default.sh\n",
        "cd "+self.CMSSW_BASE+"/src\n",
        "eval `scram runtime -sh`\n",
      ]

      # HTCondor getenv=True does not export LD_LIBRARY_PATH                                                                                                                     
      # --> added by hand in the script itself                                                                                                                                   
      if 'LD_LIBRARY_PATH' in os.environ:
         self.ExecLines += ['\n']
         self.ExecLines += ['export LD_LIBRARY_PATH='+os.environ['LD_LIBRARY_PATH']]
         self.ExecLines += ['\n\n']

      self.ConfigLines = [

        'batch_name = __BATCH_NAME__',

        'executable = __EXEC_FILE__',

        'output = __PATH__/logs/__NAME__'+'.out.'+'$(Cluster).$(Process)',
        'error  = __PATH__/logs/__NAME__'+'.err.'+'$(Cluster).$(Process)',
        'log    = __PATH__/logs/__NAME__'+'.log.'+'$(Cluster).$(Process)',

        '#arguments = ',

        'transfer_executable = True',

        'universe = vanilla',

        'getenv = True',

        'should_transfer_files   = IF_NEEDED',
        'when_to_transfer_output = ON_EXIT',

        'requirements = (OpSysAndVer == "SL6")',
        '#requirements = (OpSysAndVer == "CentOS7")',

        ' RequestMemory  =  2000',
        '+RequestRuntime = tomorrow',

        'queue',
      ]
      
    elif RunSystem == "IFIC":

      self.ExecLines = [

        "#!/bin/bash\n",
        "path=\"/lhome/ific/m/marherje/ParticleSwarmOptimization\"\n",
        "source ${path}/init_ilcsoft_v02-02-03.sh\n",
        "#source init_ilcsoft_v02-02-01.sh\n",
      ]

      # HTCondor getenv=True does not export LD_LIBRARY_PATH                                                                                            
      # --> added by hand in the script itself                                                                                                                           
      if 'LD_LIBRARY_PATH' in os.environ:
         self.ExecLines += ['\n']
         self.ExecLines += ['export LD_LIBRARY_PATH='+os.environ['LD_LIBRARY_PATH']]
         self.ExecLines += ['\n\n']

      self.ConfigLines = [

        'batch_name = __BATCH_NAME__',

        'executable = __EXEC_FILE__',

        'output = __PATH__/logs/__NAME__'+'.out.'+'$(Cluster).$(Process)',
        'error  = __PATH__/logs/__NAME__'+'.err.'+'$(Cluster).$(Process)',
        'log    = __PATH__/logs/__NAME__'+'.log.'+'$(Cluster).$(Process)',

        '#arguments = ',

        'transfer_executable = True',

        'universe = vanilla',

        'getenv = True',

        'should_transfer_files   = IF_NEEDED',
        'when_to_transfer_output = ON_EXIT',

        '#requirements = (OpSysAndVer == "SL6")',
        'requirements = (OpSysAndVer == "CentOS7")',

        '#RequestMemory  =  2000',
        '+JobFlavour = "normal"',

        'queue',
      ]

    elif RunSystem == "DESY":

      self.ExecLines = [
          
        "#!/bin/bash\n",
        "path=\"/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/analysis/ECAL_Sim/ML4PID/PSOforECALPID_4cat\"\n",
        "source ${path}/init_ilcsoft_v02-02-03.sh\n",
        "#source init_ilcsoft_v02-02-01.sh\n",
      ]
      
      if 'LD_LIBRARY_PATH' in os.environ:
        self.ExecLines += ['\n']
        self.ExecLines += ['export LD_LIBRARY_PATH='+os.environ['LD_LIBRARY_PATH']]
        self.ExecLines += ['\n\n']
        
        self.ConfigLines = [
          
        'batch_name = __BATCH_NAME__',

        'executable = __EXEC_FILE__',

        'output = __PATH__/logs/__NAME__'+'.out.'+'$(Cluster).$(Process)',
        'error  = __PATH__/logs/__NAME__'+'.err.'+'$(Cluster).$(Process)',
        'log    = __PATH__/logs/__NAME__'+'.log.'+'$(Cluster).$(Process)',

        '#arguments = ',

        'transfer_executable = True',

        'universe = vanilla',

        'getenv = True',

        'should_transfer_files   = IF_NEEDED',
        'when_to_transfer_output = ON_EXIT',

        '#requirements = (OpSysAndVer == "SL6")',
        'requirements = (OpSysAndVer == "CentOS7")',

        '#RequestMemory  =  2000',
        '#+RequestRuntime = tomorrow',

        'queue',
      ]

    else:
      print "could not set up the batch system ", self.RunSystem
      print "check QueHelper.py"
      exit(1)

    print "set up QueHelper\n"

  def GetExecLines(self):
    return self.ExecLines

  def GetRunLines(self):
    return self.RunLines

  def GetConfigLines(self):
    return self.ConfigLines

  def StartJob(self, runScript):

      jobID = None

      while jobID == None:

         ret = get_output(runScript, permissive=True)

         if ret != None:

            ret_lines = ret[0].split('\n')

            ret_lines = [_tmp for _tmp in ret_lines if _tmp != '']

            if len(ret_lines) == 2:

               jobID = ret_lines[1].split()[-1]+'0'

               if not is_float(jobID): jobID = None

         if jobID == None:

            WARNING('QueHelper.py -- StartJob: job submission failed, will wait 60sec and try again')

            time.sleep(60)

      return jobID

  def KillJob(self, jobID_str):

      ret = None

      while ret == None:

         ret = get_output('condor_rm '+jobID_str, permissive=True)

         if ret == None:

            WARNING('QueHelper.py -- KillJob: job removal failed, will wait 60sec and try again')

            time.sleep(60)

      return
