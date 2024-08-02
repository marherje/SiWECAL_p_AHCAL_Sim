#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "TRint.h"
#include "TDirectory.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TFormula.h"
#include "TTreeFormula.h"
#include "TMVA/Reader.h"
#include "TBrowser.h"

//New (remove later?)
#include "TObjString.h"

using namespace std;

TStopwatch                timer;

string                    outFileDir = "";

vector<string>            samples;
vector<string>            sampleFiles;
vector<string>            sampleTrees;
vector<string>            sampleSelections;

vector<string>            bins;
vector<string>            binTitles;
vector<string>            binSelections;

vector<string>            variables;
vector<string>            variableNames;
vector<string>            variableTypes;
vector<string>            variableSizes;

vector<string>            bdtvariables;
vector<string>            bdtvariableSizes;
vector<string>            bdtvariableWeights;
vector< vector<string> >  bdtvariableVariables;
vector< vector<string> >  bdtvariableVarNames;
vector< vector<string> >  bdtvariableVarTypes;

vector<TBranch*>          outBranches;

map<string,float*>        floatMap;
map<string,int*>          intMap;
map<string,float*>        floatBDTMap;
map<string,int*>          intBDTMap;

map<string,TMVA::Reader*> TMVAreaderMap;


void SetOutFileDir(string outdir){
  outFileDir = outdir;
}


void AddSample(string sample, string infileName, string treeName="ntp",string selection=""){
  samples.push_back(sample);
  sampleFiles.push_back(infileName);
  sampleTrees.push_back(treeName);
  sampleSelections.push_back(selection);
}


void AddBin(string name,string title, TString selection=""){
  bins.push_back(name);
  binTitles.push_back(title);
  binSelections.push_back(selection.Data());
}


void AddVar(string name, string type, string size="", string variable=""){
  if(variable != "")
    variables.push_back(variable);
  else
    variables.push_back(name);
  variableNames.push_back(name);
  variableTypes.push_back(type);
  variableSizes.push_back(size);
}


string GetBDTVarString(string* BDTVarArray, int* BDTVarIndices, int nVar){
  
  string bdtVariables = "";
  
  for(int iVar=0;iVar<nVar;iVar++){
    if(iVar!=0) bdtVariables += "&";
    bdtVariables += BDTVarArray[BDTVarIndices[iVar]];
  }
  
  return bdtVariables;
}


void AddBDTVar(string name, string size, string weights, string invariables){
  bdtvariables.push_back(name);
  bdtvariableSizes.push_back(size);
  bdtvariableWeights.push_back(weights);
  
  TObjArray* BDTvarArray = TString(invariables).Tokenize("&");
  bdtvariableVarNames.push_back(vector<string>());
  bdtvariableVarTypes.push_back(vector<string>());
  bdtvariableVariables.push_back(vector<string>());
  
  for(int iVar=0; iVar<BDTvarArray->GetEntries();iVar++){
    TObjArray* BDTvarSpecs = (((TObjString*) BDTvarArray->At(iVar))->GetString()).Tokenize(";");
    string bdtname = ((TObjString*) BDTvarSpecs->At(0))->GetString().Data();
    string bdttype = ((TObjString*) BDTvarSpecs->At(1))->GetString().Data();
    string bdtvar = "";    
    if(BDTvarSpecs->GetEntries()>2)
      bdtvar = ((TObjString*) BDTvarSpecs->At(2))->GetString().Data();
    
    bdtvariableVarNames.back().push_back(bdtname);
    bdtvariableVarTypes.back().push_back(bdttype);
    if(bdtvar!="")
      bdtvariableVariables.back().push_back(bdtvar);
    else
      bdtvariableVariables.back().push_back(bdtname);
  }
}


void InitBranch(TTree* outTree, string name, string type, string size="", int maxentries=100){
  string leaflist;
  if(size=="")
    leaflist = name+"/"+type;
  else
    leaflist = name+"["+size+"]/"+type;
  
  if(type=="F"){
    if(size=="")
      floatMap[name] = new float(0);
    else
      floatMap[name] = new float[maxentries];
    outBranches.push_back(outTree->Branch(name.c_str(),floatMap[name],leaflist.c_str()));
  }
  else if(type=="I"){
    if(size=="")
      intMap[name] = new int(0);
    else
      intMap[name] = new int[maxentries];
    outBranches.push_back(outTree->Branch(name.c_str(),intMap[name],leaflist.c_str()));
  }
  else
    cout << "unknown type " << type << endl;
}


void InitTMVABranch(TTree* outTree, string name, string size, string weights, vector<string> varnames, vector<string> types){
  TMVAreaderMap[name] = new TMVA::Reader();
  
  for(size_t iVar=0; iVar<varnames.size();iVar++){
    if(types[iVar]=="F"){
      if(floatBDTMap.count(varnames[iVar])==0)
        floatBDTMap[varnames[iVar]] = new float(0);
      TMVAreaderMap[name]->AddVariable(varnames[iVar].c_str(),floatBDTMap[varnames[iVar]]);
      cout << "BDT " << name << ": Initialized Variable " << varnames[iVar] << " (Type " << types[iVar] << ")" << endl;
    }
    else if(types[iVar]=="I"){
      if(intBDTMap.count(varnames[iVar])==0)
        intBDTMap[varnames[iVar]] = new int(0);
      TMVAreaderMap[name]->AddVariable(varnames[iVar].c_str(),intBDTMap[varnames[iVar]]);
      cout << "BDT " << name << ": Initialized Variable " << varnames[iVar] << " (Type " << types[iVar] << ")" << endl;
    }
    else
      cout << "Error! unknown variable " << varnames[iVar] << endl;
  }
  
  TMVAreaderMap[name]->BookMVA("multiclass",weights.c_str()); //Was BDT
  
  InitBranch(outTree,name,"F",size); 
}


void InitAllBranches(TTree* outTree){
  for(size_t iVar=0; iVar<variables.size();iVar++)
    InitBranch(outTree,variableNames[iVar],variableTypes[iVar],variableSizes[iVar]);
    
  for(size_t iVar=0; iVar<bdtvariables.size();iVar++)
    InitTMVABranch(outTree,bdtvariables[iVar],bdtvariableSizes[iVar],bdtvariableWeights[iVar],bdtvariableVarNames[iVar],bdtvariableVarTypes[iVar]);
}


void FillVars(TChain* inTree, int sampleIndex, int binIndex, bool appendVars = false){
  std::cout<<"filling vars"<<std::endl;
  //string outFileName = outFileDir+samples[sampleIndex]+"_"+bins[binIndex]+"_tree.root";
  string outFileName = outFileDir+"/"+samples[sampleIndex]+".root";

  TFile* outFile = new TFile(outFileName.c_str(),"recreate");
  TTree* outTree;
  
  if(appendVars){

    if(binSelections[binIndex]!="" || sampleSelections[sampleIndex]!=""){
      string selectionString= "";

      if(binSelections[binIndex]!="") selectionString+="(";
      selectionString+=binSelections[binIndex];
      if(binSelections[binIndex]!="") selectionString+=")";
      if(binSelections[binIndex]!="" && sampleSelections[sampleIndex]!="") selectionString+="&&";
      if(sampleSelections[sampleIndex]!="") selectionString+="(";
      selectionString+=sampleSelections[sampleIndex];
      if(sampleSelections[sampleIndex]!="") selectionString+=")";

      outTree = (TTree*) (inTree->CopyTree(selectionString.c_str())->Clone());
    }
    else
      outTree = inTree->CloneTree();
    //DANGERZONE
//     inTree = outTree;
  }
  else{
    
    outTree = new TTree(inTree->GetName(),inTree->GetName());
  }
  std::cout<<inTree->GetEntries()<<std::endl;
  
  TTreeFormula* sampleSelFormula = new TTreeFormula("noSampleSelection","1",inTree);
  TTreeFormula* binSelFormula = new TTreeFormula("noBinSelection","1",inTree);
  
  if(!appendVars){
    if(sampleSelections[sampleIndex]!=""){
      string sampleSelName = samples[sampleIndex]+"Selection";
      sampleSelFormula = new TTreeFormula(sampleSelName.c_str(),sampleSelections[sampleIndex].c_str(),inTree);
    }
    if(binSelections[binIndex]!=""){
      string binSelName = bins[binIndex]+"Selection";
      binSelFormula = new TTreeFormula(binSelName.c_str(),binSelections[binIndex].c_str(),inTree);
    }
  }  
  
  InitAllBranches(outTree);
  
  vector<TTreeFormula*> variableFormula;
  for(size_t iVar=0;iVar<variables.size();iVar++){
    variableFormula.push_back(new TTreeFormula(variableNames[iVar].c_str(),variables[iVar].c_str(),inTree));
    cout<<"created formula"<<variableNames[iVar].c_str()<<" "<<variables[iVar].c_str()<<" "<<inTree<<endl;
  }
  
  vector<TTreeFormula*> BDTvariableSizesFormula;
  for(size_t iVar=0;iVar<bdtvariables.size();iVar++)
    BDTvariableSizesFormula.push_back(new TTreeFormula(bdtvariableSizes[iVar].c_str(),bdtvariableSizes[iVar].c_str(),inTree->GetTree()));
  
  map<string,TTreeFormula*> BDTvariableFormula;
  for(size_t iBDTVar=0;iBDTVar<bdtvariables.size();iBDTVar++){
    for(size_t iBDTVarVar=0;iBDTVarVar<bdtvariableVarNames[iBDTVar].size();iBDTVarVar++){
      if(bdtvariableVarTypes[iBDTVar][iBDTVarVar]=="F" && floatMap.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0) continue;
      else if(bdtvariableVarTypes[iBDTVar][iBDTVarVar]=="I" && intMap.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0) continue;
      else if(BDTvariableFormula.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0) continue;
      
      BDTvariableFormula[bdtvariableVarNames[iBDTVar][iBDTVarVar]] = new TTreeFormula(bdtvariableVarNames[iBDTVar][iBDTVarVar].c_str(),bdtvariableVariables[iBDTVar][iBDTVarVar].c_str(),inTree->GetTree());
    }
  }
  
  int nEntries = inTree->GetEntries();
//   inTree->LoadTree(0);
    Int_t previous = inTree->GetTreeNumber();
  for(int iEntry=0;iEntry<nEntries;iEntry++){
    
    inTree->GetEntry(iEntry);
    Int_t current = inTree->GetTreeNumber();
    if (previous!=current)sampleSelFormula->UpdateFormulaLeaves();
    if(iEntry%10000 == 0){
      float rt = timer.RealTime();
      float cput = timer.CpuTime();
      cout << "Event " << iEntry << ": " << rt <<"s real time, " << cput << "s cpu time" << endl;
      timer.Continue();
    }
    
    if(!appendVars){
      binSelFormula->GetNdata();
      sampleSelFormula->GetNdata();
      if(!(binSelFormula->EvalInstance(0)) || !(sampleSelFormula->EvalInstance(0)))
        continue;
    }
    
    
    for(size_t iVar=0;iVar<variables.size();iVar++){
      int nInstances = variableFormula[iVar]->GetNdata();
      if (previous!=current)variableFormula[iVar]->UpdateFormulaLeaves();
      if(nInstances == 1){
        if(variableTypes[iVar]=="F"){
          *floatMap[variableNames[iVar]] = variableFormula[iVar]->EvalInstance(0);
	}
        else if(variableTypes[iVar]=="I")
          *intMap[variableNames[iVar]] = variableFormula[iVar]->EvalInstance(0);
        else
          cout << "unknown type " << variableTypes[iVar] << endl;
      }
      else{
        for(int iInstance=0;iInstance<nInstances;++iInstance){
          if(variableTypes[iVar]=="F")
            floatMap[variableNames[iVar]][iInstance] = variableFormula[iVar]->EvalInstance(iInstance);
          else if(variableTypes[iVar]=="I")
            intMap[variableNames[iVar]][iInstance] = variableFormula[iVar]->EvalInstance(iInstance);
          else
            cout << "unknown type " << variableTypes[iVar] << endl;
        }
      }      
    }
    for(size_t iBDTVar=0;iBDTVar<bdtvariables.size();iBDTVar++){
      
      BDTvariableSizesFormula[iBDTVar]->GetNdata();
      int nBTDInstances = BDTvariableSizesFormula[iBDTVar]->EvalInstance(0);
      
      for(int iInstance=0;iInstance<nBTDInstances;iInstance++){
        for(size_t iBDTVarVar=0;iBDTVarVar<bdtvariableVarNames[iBDTVar].size();iBDTVarVar++){
          if(bdtvariableVarTypes[iBDTVar][iBDTVarVar]=="F"){
            if(floatMap.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0)
              *floatBDTMap[bdtvariableVarNames[iBDTVar][iBDTVarVar]] = floatMap[bdtvariableVarNames[iBDTVar][iBDTVarVar]][iInstance];
            else if(BDTvariableFormula.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0){
              BDTvariableFormula[bdtvariableVarNames[iBDTVar][iBDTVarVar]]->GetNdata();
              *floatBDTMap[bdtvariableVarNames[iBDTVar][iBDTVarVar]] = BDTvariableFormula[bdtvariableVarNames[iBDTVar][iBDTVarVar]]->EvalInstance(iInstance);
            }
            else
              cout << "Error! Unknown Variable " << bdtvariableVariables[iBDTVar][iBDTVarVar] << endl;
          }
          else if(bdtvariableVarTypes[iBDTVar][iBDTVarVar]=="I"){
            if(intMap.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0)
              *intBDTMap[bdtvariableVarNames[iBDTVar][iBDTVarVar]] = intMap[bdtvariableVarNames[iBDTVar][iBDTVarVar]][iInstance];
            else if(BDTvariableFormula.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0){
              BDTvariableFormula[bdtvariableVarNames[iBDTVar][iBDTVarVar]]->GetNdata();
              *intBDTMap[bdtvariableVarNames[iBDTVar][iBDTVarVar]] = BDTvariableFormula[bdtvariableVarNames[iBDTVar][iBDTVarVar]]->EvalInstance(iInstance);
            }
            else
              cout << "Error! Unknown Variable " << bdtvariableVarNames[iBDTVar][iBDTVarVar] << endl;
          } 
          else
            cout << "Error! Unknown Type " << bdtvariableVarTypes[iBDTVar][iBDTVarVar] << endl;
        }
        
        floatMap[bdtvariables[iBDTVar]][iInstance] = TMVAreaderMap[bdtvariables[iBDTVar]]->EvaluateMVA("BDT");
      }
    }
    
    if(!appendVars)
      outTree->Fill();
    else
      for(size_t iBranch=0;iBranch<outBranches.size();iBranch++)
        outBranches[iBranch]->Fill();
      
    previous=current;
  }
  
  outBranches.clear();
  
  outTree->Write("",TObject::kOverwrite); //outTree->Write("", TObject::kOverwrite);
  outFile->Close();
}


void FillVarsFlat(TTree* inTree, int sampleIndex, int binIndex){
  string outFileName = outFileDir+"/"+samples[sampleIndex]+"_"+bins[binIndex]+"_tree.root";
  
  TFile* outFile = new TFile(outFileName.c_str(),"recreate");
  TTree* outTree = new TTree(inTree->GetName(),inTree->GetName());
  
  InitAllBranches(outTree);
  
  TTreeFormula* binSelFormula = new TTreeFormula("noBinSelection","1",inTree);
  if(binSelections[binIndex]!=""){
    string binSelName = bins[binIndex]+"_Selection";
    binSelFormula = new TTreeFormula(binSelName.c_str(),binSelections[binIndex].c_str(),inTree);
  }
  TTreeFormula* sampleSelFormula = new TTreeFormula("noSampleSelection","1",inTree);
  if(sampleSelections[sampleIndex]!=""){
      string sampleSelName = samples[sampleIndex]+"Selection";
      sampleSelFormula = new TTreeFormula(sampleSelName.c_str(),sampleSelections[sampleIndex].c_str(),inTree);
  }
  
  vector<TTreeFormula*> variableFormula;
  for(size_t iVar=0;iVar<variables.size();iVar++)
    variableFormula.push_back(new TTreeFormula(variableNames[iVar].c_str(),variables[iVar].c_str(),inTree));
  
  map<string,TTreeFormula*> BDTvariableFormula;
  for(size_t iBDTVar=0;iBDTVar<bdtvariables.size();iBDTVar++){
    for(size_t iBDTVarVar=0;iBDTVarVar<bdtvariableVarNames[iBDTVar].size();iBDTVarVar++){
      if(bdtvariableVarTypes[iBDTVar][iBDTVarVar]=="F" && floatMap.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0) continue;
      else if(bdtvariableVarTypes[iBDTVar][iBDTVarVar]=="I" && intMap.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0) continue;
      else if(BDTvariableFormula.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0) continue;
      
      BDTvariableFormula[bdtvariableVarNames[iBDTVar][iBDTVarVar]] = new TTreeFormula(bdtvariableVarNames[iBDTVar][iBDTVarVar].c_str(),bdtvariableVariables[iBDTVar][iBDTVarVar].c_str(),inTree);
    }
  }
  
  for(int iEntry=0;iEntry<inTree->GetEntries();iEntry++){
    inTree->GetEntry(iEntry);
    
    binSelFormula->GetNdata();
    if(!(binSelFormula->EvalInstance(0))) continue;
    
    sampleSelFormula->GetNdata();  
    for(size_t iVar=0;iVar<variables.size();iVar++)
      variableFormula[iVar]->GetNdata();
    
    for(int iInstance=0;iInstance<variableFormula[0]->GetNdata();iInstance++){
      if(!(sampleSelFormula->EvalInstance(iInstance))) continue;
      
      for(size_t iVar=0;iVar<variables.size();iVar++){
        if(variableTypes[iVar]=="F")
          *floatMap[variableNames[iVar]] = variableFormula[iVar]->EvalInstance(iInstance);
        else if(variableTypes[iVar]=="I")
          *intMap[variableNames[iVar]] = variableFormula[iVar]->EvalInstance(iInstance);
        else
          cout << "unknown type " << variableTypes[iVar] << endl;
      }
      
      for(size_t iBDTVar=0;iBDTVar<bdtvariables.size();iBDTVar++){
        for(size_t iBDTVarVar=0;iBDTVarVar<bdtvariableVarNames[iBDTVar].size();iBDTVarVar++){
          if(bdtvariableVarTypes[iBDTVar][iBDTVarVar]=="F"){
            if(floatMap.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0)
              *floatBDTMap[bdtvariableVarNames[iBDTVar][iBDTVarVar]] = *floatMap[bdtvariableVarNames[iBDTVar][iBDTVarVar]];
            else if(BDTvariableFormula.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0)
              *floatBDTMap[bdtvariableVarNames[iBDTVar][iBDTVarVar]] = BDTvariableFormula[bdtvariableVarNames[iBDTVar][iBDTVarVar]]->EvalInstance(iInstance);
            else
              cout << "Error! Unknown Variable " << bdtvariableVariables[iBDTVar][iBDTVarVar] << endl;
          }
          else if(bdtvariableVarTypes[iBDTVar][iBDTVarVar]=="I"){
            if(intMap.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0)
              *intBDTMap[bdtvariableVarNames[iBDTVar][iBDTVarVar]] = *intMap[bdtvariableVarNames[iBDTVar][iBDTVarVar]];
            else if(BDTvariableFormula.count(bdtvariableVarNames[iBDTVar][iBDTVarVar])!=0)
              *intBDTMap[bdtvariableVarNames[iBDTVar][iBDTVarVar]] = BDTvariableFormula[bdtvariableVarNames[iBDTVar][iBDTVarVar]]->EvalInstance(iInstance);
            else
              cout << "Error! Unknown Variable " << bdtvariableVarNames[iBDTVar][iBDTVarVar] << endl;
          } 
          else
            cout << "Error! Unknown Type " << bdtvariableVarTypes[iBDTVar][iBDTVarVar] << endl;  
        }
          
        *floatMap[bdtvariables[iBDTVar]] = TMVAreaderMap[bdtvariables[iBDTVar]]->EvaluateMVA("multiclass"); //Was "BDT"
      }
      
      outTree->Fill();
    }
  }
  
  outTree->Write();
  outFile->Close();
}


void FillAllTrees(bool flatTrees = false, bool appendVars = false){
  timer.Start();
  
  for(size_t iSample=0;iSample<samples.size();iSample++){
    
//     TFile* inFile = new TFile(sampleFiles[iSample].c_str());
//     TTree* inTree = (TTree*) inFile->Get(sampleTrees[iSample].c_str());
    std::cout<<"adding Trees to chain"<<std::endl;
    TChain* inTree = new TChain(sampleTrees[iSample].c_str());
    inTree->Add(sampleFiles[iSample].c_str());
    
    for( size_t iBin=0; iBin<bins.size();iBin++){
      cout << "FillAllTrees Start Sample " << samples[iSample] << " and Bin " << bins[iBin] << endl;
      if(!flatTrees){
        FillVars(inTree,iSample,iBin,appendVars);
      }
      else{
        FillVarsFlat(inTree,iSample,iBin);
      }
    }

//     inFile->Close();
  }
  
  timer.Stop();
  float rt = timer.RealTime();
  float cput = timer.CpuTime();
  cout << "Tree Skimming took " << rt <<"s real time, " << cput << "s cpu time" << endl;
}
