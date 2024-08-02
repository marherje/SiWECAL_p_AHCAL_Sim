#include "TreeHeader.h"
#include <fstream>

void PrepareTrees(string configName="", const std::string& OutPath="InitData"){

string BFile;
string BTree;
string CFile;
string CTree;
string OFile;
string OTree;
string Selection="";

bool   UseFixedTrainTestSplitting=false;
string UseFixedTrainTestSplitting_Train;

// read config file
std::ifstream config(configName.c_str());
bool readline=true;
std::string linebuffer;
std::string searchedString;

while(readline){
/*  config>>linebuffer;*/
  std::getline(config,linebuffer);
  if(config.eof())readline=false;
  if(linebuffer.find("#")!=std::string::npos)continue;
//   if(linebuffer=="")std::cout<<"empty"<<std::endl;
  
  searchedString="SourceBFile=";
  if(linebuffer.find(searchedString)!=std::string::npos){
    BFile=linebuffer.substr(searchedString.length());
    std::cout<<"found BFile"<<std::endl;
  }
  searchedString="SourceBTree=";
  if(linebuffer.find(searchedString)!=std::string::npos){
    BTree=linebuffer.substr(searchedString.length());
  }
  searchedString="SourceCFile=";
  if(linebuffer.find(searchedString)!=std::string::npos){
    CFile=linebuffer.substr(searchedString.length());
    std::cout<<"found CFile"<<std::endl;
  }
  searchedString="SourceCTree=";
  if(linebuffer.find(searchedString)!=std::string::npos){
    CTree=linebuffer.substr(searchedString.length());
  }
  searchedString="SourceOFile=";
  if(linebuffer.find(searchedString)!=std::string::npos){
    OFile=linebuffer.substr(searchedString.length());
    std::cout<<"found OFile"<<std::endl;
  }
  searchedString="SourceOTree=";
  if(linebuffer.find(searchedString)!=std::string::npos){
    OTree=linebuffer.substr(searchedString.length());
  }
  searchedString="SelectionString=";
  if(linebuffer.find(searchedString)!=std::string::npos){
    Selection=linebuffer.substr(searchedString.length());
  }
  searchedString="UseFixedTrainTestSplitting=1";
  if(linebuffer.find(searchedString) != std::string::npos)
  {
    UseFixedTrainTestSplitting = true;
  }
  searchedString="UseFixedTrainTestSplitting_Train=";
  if(linebuffer.find(searchedString) != std::string::npos)
  {
    UseFixedTrainTestSplitting_Train = linebuffer.substr(searchedString.length());

    if(UseFixedTrainTestSplitting_Train.compare(0, 1, "\"") == 0)
    {
      UseFixedTrainTestSplitting_Train = UseFixedTrainTestSplitting_Train.substr(1);
    }

    if(UseFixedTrainTestSplitting_Train.compare(UseFixedTrainTestSplitting_Train.size()-1, 1, "\"") == 0)
    {
      UseFixedTrainTestSplitting_Train = UseFixedTrainTestSplitting_Train.substr(0, UseFixedTrainTestSplitting_Train.size()-1);
    }
  }

  //look for InitialVariables
  searchedString="InitialVariables";
  if(linebuffer.find(searchedString)!=std::string::npos){
    bool readingVars=true;
//     std::cout<<"InitialVariables:"<<std::endl;
    while(readingVars){
      std::getline(config,linebuffer);
      if(linebuffer.find("#")!=std::string::npos || linebuffer=="")continue;
      if(linebuffer.find("EndVariables")!=std::string::npos || config.eof()){
        readingVars=false; 
        continue;
        }
      string name,type,size,expression;
      size_t id1, id2;
      id1=linebuffer.find("\"");
      id2=linebuffer.find("\"",id1+1);
      name=linebuffer.substr(id1+1,id2-id1-1);
      id1=linebuffer.find("\"",id2+1);
      id2=linebuffer.find("\"",id1+1);
      type=linebuffer.substr(id1+1,id2-id1-1);
      id1=linebuffer.find("\"",id2+1);
      id2=linebuffer.find("\"",id1+1);
      size=linebuffer.substr(id1+1,id2-id1-1);
      id1=linebuffer.find("\"",id2+1);
      id2=linebuffer.find("\"",id1+1);
      expression=linebuffer.substr(id1+1,id2-id1-1);
      std::cout<<"Did Variable "<<name<<" "<<type<<" "<<size<<" "<<expression<<std::endl;
//       if(size=="")std::cout<<"zero"<<std::endl;
      AddVar(name,type,size,expression );
    }
  }
        
  searchedString="AdditionalVariables";
  if(linebuffer.find(searchedString)!=std::string::npos){
    bool readingVars=true;
//     std::cout<<"AdditionalVariables:"<<std::endl;
    while(readingVars){
      std::getline(config,linebuffer);
      if(linebuffer.find("#")!=std::string::npos || linebuffer=="")continue;
      if(linebuffer.find("EndVariables")!=std::string::npos || config.eof()){
        readingVars=false; 
        continue;
        }
      string name,type,size,expression;
      size_t id1, id2;
      id1=linebuffer.find("\"");
      id2=linebuffer.find("\"",id1+1);
      name=linebuffer.substr(id1+1,id2-id1-1);
      id1=linebuffer.find("\"",id2+1);
      id2=linebuffer.find("\"",id1+1);
      type=linebuffer.substr(id1+1,id2-id1-1);
      id1=linebuffer.find("\"",id2+1);
      id2=linebuffer.find("\"",id1+1);
      size=linebuffer.substr(id1+1,id2-id1-1);
      id1=linebuffer.find("\"",id2+1);
      id2=linebuffer.find("\"",id1+1);
      expression=linebuffer.substr(id1+1,id2-id1-1);
      std::cout<<"Did Variable "<<name<<" "<<type<<" "<<size<<" "<<expression<<std::endl;
//       if(size=="")std::cout<<"zero"<<std::endl;
      AddVar(name,type,size,expression );
    }
  }

  searchedString="WeightVariables";
  if(linebuffer.find(searchedString)!=std::string::npos){
    bool readingVars=true;
//     std::cout<<"WeightVariables:"<<std::endl;
    while(readingVars){
      std::getline(config,linebuffer);
      if(linebuffer.find("#")!=std::string::npos || linebuffer=="")continue;
      if(linebuffer.find("EndVariables")!=std::string::npos || config.eof()){
        readingVars=false; 
        continue;
        }
      string name,type,size,expression;
      size_t id1, id2;
      id1=linebuffer.find("\"");
      id2=linebuffer.find("\"",id1+1);
      name=linebuffer.substr(id1+1,id2-id1-1);
      id1=linebuffer.find("\"",id2+1);
      id2=linebuffer.find("\"",id1+1);
      type=linebuffer.substr(id1+1,id2-id1-1);
      id1=linebuffer.find("\"",id2+1);
      id2=linebuffer.find("\"",id1+1);
      size=linebuffer.substr(id1+1,id2-id1-1);
      id1=linebuffer.find("\"",id2+1);
      id2=linebuffer.find("\"",id1+1);
      expression=linebuffer.substr(id1+1,id2-id1-1);
      std::cout<<"Did Variable "<<name<<" "<<type<<" "<<size<<" "<<expression<<std::endl;
//       if(size=="")std::cout<<"zero"<<std::endl;
      AddVar(name,type,size,expression );
    }
  }

  searchedString="SpectatorVariables";
  if(linebuffer.find(searchedString)!=std::string::npos){
    bool readingVars=true;
//     std::cout<<"WeightVariables:"<<std::endl;
    while(readingVars){
      std::getline(config,linebuffer);
      if(linebuffer.find("#")!=std::string::npos || linebuffer=="")continue;
      if(linebuffer.find("EndVariables")!=std::string::npos || config.eof()){
        readingVars=false; 
        continue;
        }
      string name,type,size,expression;
      size_t id1, id2;
      id1=linebuffer.find("\"");
      id2=linebuffer.find("\"",id1+1);
      name=linebuffer.substr(id1+1,id2-id1-1);
      id1=linebuffer.find("\"",id2+1);
      id2=linebuffer.find("\"",id1+1);
      type=linebuffer.substr(id1+1,id2-id1-1);
      id1=linebuffer.find("\"",id2+1);
      id2=linebuffer.find("\"",id1+1);
      size=linebuffer.substr(id1+1,id2-id1-1);
      id1=linebuffer.find("\"",id2+1);
      id2=linebuffer.find("\"",id1+1);
      expression=linebuffer.substr(id1+1,id2-id1-1);
      std::cout<<"Did Variable "<<name<<" "<<type<<" "<<size<<" "<<expression<<std::endl;
//       if(size=="")std::cout<<"zero"<<std::endl;
      AddVar(name,type,size,expression );
    }
  }

}
config.close();

if(UseFixedTrainTestSplitting)
{
  string selTrain = " ("+UseFixedTrainTestSplitting_Train+")";
  string selTest  = "!("+UseFixedTrainTestSplitting_Train+")";

  if(Selection != "")
  {
    selTrain += " * ("+Selection+")";
    selTest  += " * ("+Selection+")";
  }

  std::cout << " Train = \"" << selTrain << "\"" << std::endl;

  AddSample("B_Train", BFile, BTree, selTrain);
  AddSample("B_Test" , BFile, BTree, selTest );

  AddSample("C_Train", CFile, CTree, selTrain);
  AddSample("C_Test" , CFile, CTree, selTest );

  AddSample("O_Train", OFile, OTree, selTrain);
  AddSample("O_Test" , OFile, OTree, selTest );
}
else
{
  AddSample("B",BFile,BTree,Selection);
  AddSample("C",CFile,CTree,Selection);
  AddSample("O",OFile,OTree,Selection);
}

std::cout<<"added samples"<<std::endl;

 AddBin("all","","");

bool flatTrees = false;
bool appendVariables = false;

SetOutFileDir(OutPath);
FillAllTrees(flatTrees,appendVariables);

}
