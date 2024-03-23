import sys
from subprocess import call
from glob import glob
import os

PSOdirs=[
	]

PSOdirs=sys.argv[1:]

listOfPSOs=[]
listOfPrelPSOs=[]

for d in PSOdirs:
  for p in glob(d+"/*"):
    if os.path.isdir(p):
      if os.path.isfile(p+"/PSOResult.txt"):
        listOfPrelPSOs.append([d,p.rsplit("/",1)[1]])
headerrow=[["path","PSO","Iteration","ProblemFrac","bestROC","bestKS"]]
for pso in listOfPrelPSOs:
  listOfOutFilesP0=glob(pso[0]+"/"+pso[1]+"/Particles/Particle0/logs/*.o*")
  listOfErrorFiles=glob(pso[0]+"/"+pso[1]+"/Particles/Particle*/logs/*.e*")
  nIts=len(listOfOutFilesP0)
  nErrors=0
  for ef in listOfErrorFiles:
    if os.stat(ef).st_size>0:
      nErrors+=1
  bestROC=0
  bestKS=0
  resfile=open(pso[0]+"/"+pso[1]+"/PSOResult.txt","r")
  reslist=list(resfile)
  bestline=reslist[9]
  bestROC=bestline.split(",")[0]
  bestKS=bestline.split(",")[1]
  resfile.close()
  bestROC=bestROC.replace("[","")
  thispso=[]
  thispso.append(pso[0])
  thispso.append(pso[1])
  thispso.append(str(nIts))
  thispso.append(str(nErrors)+"/"+str(len(listOfErrorFiles)))
  thispso.append(str(bestROC))
  thispso.append(str(bestKS))
  listOfPSOs.append(thispso)
  
assert(len(headerrow[0])==len(listOfPSOs[0]))
maxColLengths=[0]*len(listOfPSOs[0])

for pso in headerrow+listOfPSOs:
  for ic, c in enumerate(pso):
    if len(c)>maxColLengths[ic]:
      maxColLengths[ic]=len(c)

for pso in headerrow+listOfPSOs:
  string=""
  for ic, c in enumerate(pso):
    ss=c
    while len(ss)<maxColLengths[ic]:
      ss+=" "
    string+=ss
    string+="\t"
  string+="\n"
  print string
  
  
  
