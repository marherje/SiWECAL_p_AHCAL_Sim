#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>

#include "TGraph.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TMath.h"

#include "TStopwatch.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
//include "TMVA/MethodBDT.h"
#include "TMVA/ResultsClassification.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TVirtualFitter.h"

#include "TMVA/DataSet.h"
#include "TMVA/ResultsMulticlass.h"
#include "TMatrixT.h"
//include "TMVA/GetMulticlassConfusionMatrix.h"


Double_t TimeTotal=0.0;
Double_t TimeTreePreparation=0.0;
Double_t TimeTraining=0.0;
Double_t TimeTesting=0.0;
int NTrainings=3;




class BDTVar :public TObject {
  public:
    BDTVar();
    ~BDTVar();
    
    TString name;
//     Double_t FOMI;
//     Double_t KS;
//     Double_t Separation;
    
//     ClassDef(BDTVar,1);
    
};

BDTVar::BDTVar(){
}

BDTVar::~BDTVar(){
}


void DoTraining(std::vector<BDTVar*> UsedVars,std::vector<BDTVar*> UnUsedVars,TString FOMType, TString FactoryString, TString PrepString, TString BWeight, TString CWeight, TString OWeight, TString BTreeName, TString CTreeName,TString OTreeName ,TString MethodType, TString MethodString, int particleNumber, Double_t* testFOM, Double_t* testKS,Int_t UseFixedTrainTestSplitting, TString PlotName="NONE" ){
   std::cout<<"----------------------------------------------------------------"<<std::endl;

   //prepare TMVA
   TMVA::IMethod* im;
   //timing studies

   TStopwatch* thisTimer = new TStopwatch();

   Double_t thisTreeTime = 0.0;
   Double_t thisTrainingTime = 0.0;
   Double_t thisTestingTime = 0.0;

   TMVA::Tools::Instance();

   TString outfileName("TMVAMulticlass.root");

   int nUsedVars=UsedVars.size();

   TFile* outputFile = TFile::Open(outfileName, "RECREATE");

   TMVA::Factory *factory = new TMVA::Factory("TMVAMulticlass", outputFile, FactoryString); //TMVAMulticlass
   std::cout << "MethodString = " << MethodString << std::endl;

   TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataloader");

   std::cout<<"Using Variables"<<std::endl;
   for(int k=0; k<UsedVars.size(); k++)
   {
      if (UsedVars.at(k)->name.Contains("Combination"))
      {
        std::cout<<"Adding Spectator " <<k<<" "<<UsedVars.at(k)->name<<std::endl; 
        dataloader->AddSpectator(UsedVars.at(k)->name);
      }
      else
      {
        std::cout<<k<<" "<<UsedVars.at(k)->name<<std::endl; 
        dataloader->AddVariable(UsedVars.at(k)->name);
      }
     
   }
   //   dataloader->AddSpectator("nvtx","","I");
   //   dataloader->AddSpectator("nvtxall","","I");
   //   dataloader->AddSpectator( "correctCombination", "", "", 0, 1, "I" );
   //   dataloader->AddSpectator( "swappedCombination", "", "", 0, 1, "I" );
   thisTimer->Start();

   TFile *inputBTrain;
   TFile *inputCTrain;
   TFile *inputOTrain;
   TTree *BTrain;
   TTree *CTrain;
   TTree *OTrain;
   TFile *inputBTest;
   TFile *inputCTest;
   TFile *inputOTest;
   TTree *BTest;
   TTree *CTest;
   TTree *OTest;
   TFile *inputB;
   TFile *inputC;
   TFile *inputO;
   TTree *TreeB;
   TTree *TreeC;
   TTree *TreeO;

   Double_t wBfChi2=0.0;
   Double_t wCfChi2=0.0;
   Double_t wOfChi2=0.0;

   //TCut CutNaN="trk1d0sig!=0&&!TMath::IsNaN(trk1z0sig)&&!TMath::IsNaN(trk2z0sig)&&!TMath::IsNaN(jprobz25sigma)&&!TMath::IsNaN(jprobr25sigma)&&!TMath::IsNaN(jprobr2)&&!TMath:IsNaN(jprobz2)";

   if(UseFixedTrainTestSplitting == 1)
   {
     inputBTrain = TFile::Open("B_Train.root");
     inputCTrain = TFile::Open("C_Train.root");
     inputOTrain = TFile::Open("O_Train.root");

     BTrain = (TTree*) inputBTrain->Get(BTreeName);
     CTrain = (TTree*) inputCTrain->Get(CTreeName);
     OTrain = (TTree*) inputOTrain->Get(OTreeName);

     
     dataloader->AddTree(BTrain,"B_Train",1.0,"",TMVA::Types::kTraining);
     dataloader->AddTree(CTrain,"C_Train", 1.0,"",TMVA::Types::kTraining);                             
     dataloader->AddTree(OTrain,"O_Train", 1.0,"",TMVA::Types::kTraining);

     inputBTest = TFile::Open("B_Test.root");
     inputCTest = TFile::Open("C_Test.root");
     inputOTest = TFile::Open("O_Test.root");

     BTest = (TTree*) inputBTest->Get(BTreeName);
     CTest = (TTree*) inputCTest->Get(CTreeName);
     OTest = (TTree*) inputOTest->Get(OTreeName);

     dataloader->AddTree(BTest,"B_Test",1.0,"",TMVA::Types::kTesting); 
     dataloader->AddTree(CTest,"C_Test",1.0,"",TMVA::Types::kTesting);
     dataloader->AddTree(OTest,"O_Test",1.0,"",TMVA::Types::kTesting);

     wBfChi2 = BTest->GetMinimum(CWeight) * BTest->GetEntries();
     wCfChi2 = CTest->GetMinimum(CWeight) * CTest->GetEntries();
     wOfChi2 = OTest->GetMinimum(OWeight) * OTest->GetEntries();
   }
   else
   {
     inputB = TFile::Open( "B.root" );
     inputC = TFile::Open( "C.root" );
     inputO = TFile::Open( "O.root" );

     TreeB = (TTree*) inputB->Get(BTreeName);
     TreeC = (TTree*) inputC->Get(CTreeName);
     TreeO = (TTree*) inputO->Get(OTreeName);

     wBfChi2 = (TreeB->GetMinimum(BWeight) * TreeB ->GetEntries());
     wCfChi2 = (TreeC->GetMinimum(CWeight) * TreeC ->GetEntries());
     wOfChi2 = (TreeO->GetMinimum(OWeight) * TreeO ->GetEntries());

     dataloader->AddTree(TreeB,"jet_B",1.0,"");
     dataloader->AddTree(TreeC,"jet_C" ,1.0,"");
     dataloader->AddTree(TreeO,"jet_O" ,1.0,"");
   }

   dataloader->SetWeightExpression( BWeight,"jet_B");
   dataloader->SetWeightExpression( CWeight,"jet_C");
   dataloader->SetWeightExpression( OWeight,"jet_O");

   std::cout << "Running TMVA::DataLoader::PrepareTrainingAndTestTree(\"\", \"\", \""+PrepString+"\")" << std::endl;
   dataloader->PrepareTrainingAndTestTree("", PrepString);
      
   //check method and book it
   if(MethodType=="TMVA::Types::kBDT")
   {
     factory->BookMethod(dataloader, TMVA::Types::kBDT,"myMVA", MethodString);
   }
   else if(MethodType=="TMVA::Types::kLikelihood")
   {
     factory->BookMethod(dataloader, TMVA::Types::kLikelihood,"myMVA", MethodString);
   }
   else if(MethodType=="TMVA::Types::kMLP")
   {
     factory->BookMethod(dataloader, TMVA::Types::kMLP,"myMVA", MethodString);
   }
   else
   {
     std::cout<<"dont know "<<MethodType<<std::endl;
     exit(1);
   }

   //timing studies
   thisTreeTime=thisTimer->RealTime();
   thisTimer->Start();

   std::cout << "Running TMVA::Factory::TrainAllMethods()" << std::endl;
   factory->TrainAllMethods();

   std::cout<<"--------Training Done------"<<std::endl;

   //timing studies
   thisTrainingTime=thisTimer->RealTime();
   thisTimer->Start();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();
   factory->EvaluateAllMethods();


   // Weight - N. of events for each class, to use for AUC:
   float Total_B=TreeB->GetEntries();
   float Total_C=TreeC->GetEntries();
   float Total_O=TreeO->GetEntries();
   float Total_All=Total_B+Total_C+Total_O;
   
   float Weight_B=Total_B/Total_All;
   float Weight_C=Total_C/Total_All;
   float Weight_O=Total_O/Total_All;


   double ROC_B=factory->GetROCIntegral(dataloader,"myMVA",0); //MethodString? or "myMVA"?
   double ROC_C=factory->GetROCIntegral(dataloader,"myMVA",1);
   double ROC_O=factory->GetROCIntegral(dataloader,"myMVA",2);

   std::cout<<"ROC Integrals(B,C,O): "<<ROC_B<<" | "<<ROC_C<<" | "<<ROC_O<<std::endl;
   
   // COMENTAR DESDE AQUI
   double Weighted_AUC=(Weight_B*ROC_B+Weight_C*ROC_C+Weight_O*ROC_O);
   std::cout<<"FOM (Weighted ROC Integral sum): "<<Weighted_AUC<<std::endl;
   
   double Weighted_1_AUC=(ROC_B+ROC_C+ROC_O)/3;
   std::cout<<"FOM (Weighted over 1 ROC Integral): "<<Weighted_1_AUC<<std::endl;
   

   //im = factory->GetMethod(dataloader->GetName(), "myMVA");
   //im = factory->GetMethod(MethodType,"myMVA");
   //TMVA::MethodBase* myMethod = dynamic_cast<TMVA::MethodBase*>(im);
   //Double_t err;

   //TMatrixT<double> ConffusionMatrixTry=myMethod->TMVA::MethodBase::GetMulticlassConfusionMatrix(0.3,TMVA::Types::kTraining);
   //dumpElements(ConffusionMatrixTry);
   
   //UInt_t TryEventsDeleteLater=myMethod->GetNEvents();
   //std::cout<<"Reading entries from Method: "<<TryEventsDeleteLater<<std::endl;
   //myMethod->TestClassification();
   
   //Double_t ROCC1 = myMethod->GetROCIntegral();
   //Double_t ROC = myMethod->GetEfficiency();
   
   //std::cout<<"ROCC1: "<<ROCC1<<std::endl;
   
   //buff is needed for KS test
   /*   Double_t buff = myMethod->GetTrainingEfficiency("Efficiency:0.5");

   std::cout << "doing sep" << std::endl;

   Double_t sep = myMethod->GetSeparation();*/
//   Double_t sep=0.0;

   /*TMVA::DataSet* myData = myMethod->Data();*/
   //myData->SetCurrentType(TMVA::Types::kTesting);
   
   /*
   TMVA::ResultsClassification* mvaRes = dynamic_cast<TMVA::ResultsClassification*>(myData->GetResults(myMethod->GetMethodName()));
   
   TH1* mva_s = (TH1*) mvaRes->GetHist("MVA_S");
   TH1* mva_b = (TH1*) mvaRes->GetHist("MVA_B");
   TH1* mva_effS = (TH1*) mvaRes->GetHist("MVA_EFF_S");
   TH1* mva_effB = (TH1*) mvaRes->GetHist("MVA_EFF_B");

   TH1D* mva_s_high = (TH1D*) mvaRes->GetHist("MVA_HIGHBIN_S");
   TH1D* mva_b_high = (TH1D*) mvaRes->GetHist("MVA_HIGHBIN_B");
   */
   // DESCOMENTAR HASTA AQUI
/*   TH1* mva_FOMCurve = (TH1*) mvaRes->GetHist("MVA_myMVA_rejBvsS");*/

   // DESCOMENTAR DESDE AQUI
   /*
   //do chi2 stuff
   Double_t chi2FOM=0.0;
   if(FOMType=="Chi2_B_muSB"){
     chi2FOM=GetChi2FOM(mva_s_high,wSfChi2,mva_b_high,wBfChi2,"Chi2_B_muSB");
   }

   int roccNBins=mva_effS->GetNbinsX();
   TGraph* mva_ROCCurve_H = new TGraph(roccNBins);
   for(int ib=1;ib<=roccNBins;ib++){
     Double_t eS=mva_effS->GetBinContent(ib);
     Double_t rB=1.0-mva_effB->GetBinContent(ib);
     mva_ROCCurve_H->SetPoint(ib,eS,rB);
   }
   mva_ROCCurve_H->Sort();
   if(PlotName!="NONE"){
     mva_ROCCurve_H->SaveAs(PlotName+".root");
   }

   int bin=0;
   bin=mva_effS->FindLastBinAbove(0.1)+1;
   Double_t effS_01 = mva_effS->GetBinContent(bin);
   Double_t rejB_01 = 1.0-mva_effB->GetBinContent(bin);
//    std::cout<<bin<<std::endl;
   std::cout<<"rejB at 0.1 effS "<<rejB_01<<" at "<<effS_01<<std::endl;
   bin=mva_effS->FindLastBinAbove(0.3)+1;
//    std::cout<<bin<<std::endl;
   Double_t effS_03 = mva_effS->GetBinContent(bin);
   Double_t rejB_03 = 1.0-mva_effB->GetBinContent(bin);
   std::cout<<"rejB at 0.3 effS "<<rejB_03<<" at "<<effS_03<<std::endl;
   bin=mva_effS->FindLastBinAbove(0.5)+1;
//    std::cout<<bin<<std::endl;
   Double_t effS_05 = mva_effS->GetBinContent(bin);
   Double_t rejB_05 = 1.0-mva_effB->GetBinContent(bin);
   std::cout<<"rejB at 0.5 effS "<<rejB_05<<" at "<<effS_05<<std::endl;
   */

   // TRYING TO MAKE THE KS TEST:
   // NEW PART:
   //im = factory->GetMethod(dataloader->GetName(), "myMVA");
   /*
   std::cout<<"DEBUG for KS"<<std::endl;
   im = factory->GetMethod(MethodType,"myMVA");  
   std::cout<<"passed im"<<std::endl;
   TMVA::MethodBase* myMethod = dynamic_cast<TMVA::MethodBase*>(im);
   std::cout<<"passed MethodBase"<<std::endl;
   TMVA::DataSet* myData = myMethod->Data();
   myData->SetCurrentType(TMVA::Types::kTesting);
   std::cout<<"passed DataSet"<<std::endl;
   
   TMVA::ResultsMulticlass* mvaRes = dynamic_cast<TMVA::ResultsMulticlass*>(myData->GetResults(myMethod->GetMethodName(),TMVA::Types::kTesting, TMVA::Types::kMulticlass));
   std::cout<<"passed ResultsClassification"<<std::endl;
   
   TH1* mva_b = (TH1*) mvaRes->GetHist("MVA_B");          
   TH1* mva_c = (TH1*) mvaRes->GetHist("MVA_C");
   TH1* mva_q = (TH1*) mvaRes->GetHist("MVA_O");
   std::cout<<"We've got the histos!"<<std::endl;

//calculate KS score
   Double_t ks_b=0;
   Double_t ks_c=0;
   Double_t ks_q=0;
   Double_t KS_BC=0;
   Double_t KS_BQ=0;
   Double_t KS_CQ=0;

   Double_t KS=0;
   //myData->SetCurrentType(TMVA::Types::kTraining);

   //TMVA::ResultsClassification* mvaResTrain = dynamic_cast<TMVA::ResultsClassification*>(myData->GetResults(myMethod->GetMethodName(),TMVA::Types::kTraining, TMVA::Types::kClassification));

   TMVA::ResultsMulticlass* mvaResTrain = dynamic_cast<TMVA::ResultsMulticlass*>(myData->GetResults(myMethod->GetMethodName(),TMVA::Types::kTraining, TMVA::Types::kMulticlass));
   TH1* mva_b_tr = (TH1*) mvaRes->GetHist("MVA_TRAIN_B");
   TH1* mva_c_tr = (TH1*) mvaRes->GetHist("MVA_TRAIN_C");
   TH1* mva_q_tr = (TH1*) mvaRes->GetHist("MVA_TRAIN_O");
   ks_b=mva_b->KolmogorovTest(mva_b_tr, "X");
   ks_c=mva_c->KolmogorovTest(mva_c_tr, "X");
   ks_q=mva_q->KolmogorovTest(mva_q_tr, "X");
   KS_BC=TMath::Min(ks_b, ks_c);
   KS_BQ=TMath::Min(ks_b, ks_q);
   KS_CQ=TMath::Min(ks_c, ks_q);
 
   Double_t auxmin1=TMath::Min(KS_BC,KS_BQ);
   Double_t auxmin2=TMath::Min(KS_BQ,KS_CQ);

   KS = TMath::Min(auxmin1,auxmin2);
   
   std::cout<<"KS b-c: " << KS_BC << std::endl;
   std::cout<<"KS b-q: " << KS_BQ << std::endl;
   std::cout<<"KS c-q: " << KS_CQ << std::endl;
   std::cout<<"KS: "     << KS  << std::endl;
   //std::cout<<"ROC: "    << ROC << std::endl;
   //std::cout<<"Sep: "    << sep << std::endl;
   */
   
   Double_t KS=0.1;
    *testKS=KS;
   //choose wanted FOM
   std::cout<<"chosen FOM "<<FOMType<<std::endl;
   if(FOMType=="Weighted_AUC"){
     *testFOM=Weighted_AUC;
   }
   if(FOMType=="Weighted_1_AUC"){
     *testFOM=Weighted_1_AUC;
   }
   else{
     std::cout<<"ERROR: Dont know the requested FOM"<<std::endl<<"Implement it !"<<std::endl;
     exit(1);
   }
   std::cout<<"KS "<<*testKS<<std::endl;
   std::cout<<"FOM "<<*testFOM<<std::endl;
   

   //timing studies
   thisTestingTime=thisTimer->RealTime();
   std::cout<<"*******************************************************"<<std::endl;
   std::cout<<"* this Iteration Tree processing time "<<thisTreeTime<<std::endl;
   std::cout<<"* this Iteration Training time "<<thisTrainingTime<<std::endl;
   std::cout<<"* this Iteration Testing  time "<<thisTestingTime<<std::endl;
   std::cout<<"*******************************************************"<<std::endl;
   TimeTreePreparation+=thisTreeTime;
   TimeTraining+=thisTrainingTime;
   TimeTesting+=thisTestingTime;

   delete thisTimer;

// free all the memory
// hopefully nothing is missed

//    delete mvaRes;
//    delete mvaResTrain;

//!!   delete myData;

//    delete mva_s;
//    delete mva_b;

//    delete mva_effS;
//    delete mva_effB;

   //delete mva_ROCCurve_H;

//    delete mva_s_tr;
//    delete mva_b_tr;

   delete factory;
   delete dataloader;

   if(UseFixedTrainTestSplitting == 0)
   {
//!!     delete B;
//!!     delete C;
//!!     delete O;

     inputB->Close();
     inputC->Close();
     inputO->Close();

//!!     delete inputS;
//!!     delete inputB;
   }
   else
   {
//!!     delete BTrain;
//!!     delete CTrain;
//!!     delete CTrain;

     inputBTrain->Close();
     inputCTrain->Close();
     inputOTrain->Close();

//!!     delete inputSTrain;
//!!     delete inputBTrain;

//!!     delete BTest;
//!!     delete CTest;

     inputBTest->Close();
     inputCTest->Close();
     inputOTest->Close();
//!!     delete inputSTest;
//!!     delete inputBTest;
   }

   outputFile->Close();
//!!   delete outputFile;
}

void Particle()
{
  TStopwatch* timerTotal=new TStopwatch();
  timerTotal->Start();

 //set up variables
 std::vector<BDTVar*> InitialVars;
 std::vector<BDTVar*> OtherVars;
 int particleNumber;
 TString Iteration="0";
 TString FOMType="Weighted_1_AUC";
 TString SaveTrainingsToTrees="True";
 Double_t KSThreshold=0.0;
 TString FactoryString="";
 TString PrepString="";
 TString BWeight="";
 TString CWeight="";
 TString OWeight="";
 TString BTreeName="ntp";
 TString CTreeName="ntp";
 TString OTreeName="ntp";
 Int_t FindBestVariables=1;
 Int_t MaxVariablesInCombination=10;
 Int_t MinVariablesInCombination=6;
 Double_t ImprovementThreshold=1.0;
 int RepeatTrainingNTimes=0;
 TString MethodType="";
 TString MethodString="";
 std::vector<TString> coordNames;
 std::vector<Double_t> coordValues;
 Double_t dumpVal;
 TString dumpName;
 Int_t UseFixedTrainTestSplitting = 0;

  //read Config File
  std::ifstream config("ParticleConfig.txt");
  TString dump="";
  
  int count=0;
  bool readline=true;
  std::cout<<"reading Config"<<std::endl;
  while(readline==true && count<400){
    count++;
    config>>dump;
    if(config.eof())readline=false;
    if(dump=="particleNumber"){
      config>>particleNumber;}
    if(dump=="Iteration"){
      config>>Iteration;}
    if(dump=="FOM"){
      config>>FOMType;}
    if(dump=="SaveTrainingsToTrees"){
      config>>SaveTrainingsToTrees;}
    if(dump=="KSThreshold"){
      config>>KSThreshold;}
    if(dump=="FactoryString"){
      config>>FactoryString;}
    if(dump=="PreparationString"){
      config>>PrepString;}
    if(dump=="BWeightExpression"){
      config>>BWeight;}
    if(dump=="CWeightExpression"){
      config>>CWeight;}
    if(dump=="OWeightExpression"){
      config>>OWeight;}
    if(dump=="BTreeName"){
      config>>BTreeName;}
    if(dump=="CTreeName"){
      config>>CTreeName;}
    if(dump=="OTreeName"){
      config>>OTreeName;}
    if(dump=="UseFixedTrainTestSplitting"){
      config>>UseFixedTrainTestSplitting;}
    if(dump=="FindBestVariables"){
      config>>FindBestVariables;}
    if(dump=="MaxVariablesInCombination"){
      config>>MaxVariablesInCombination;}
    if(dump=="ImprovementThreshold"){
      config>>ImprovementThreshold;}
    if(dump=="RepeatTrainingNTimes"){
      config>>RepeatTrainingNTimes;}
    if(dump=="MethodType"){
      config>>MethodType;}
    if(dump=="MethodParameters"){
      config>>MethodString;}
    if(dump=="coord"){
      config>>dumpName;
      config>>dumpVal;
      std::cout<<dump<<" "<<dumpName<<" "<<dumpVal<<std::endl;
      coordNames.push_back(dumpName);
      coordValues.push_back(dumpVal);}
    if(dump=="--InitialVariables--"){
      bool readVars=true;
      while(readVars==true){
        config>>dump;
        if(dump=="--EndInitVars--"){ readVars=false; continue;}
        InitialVars.push_back(new BDTVar);
        InitialVars.back()->name=dump;
      }
    }
    if(dump=="--AdditionalVariables--"){
      bool readVars=true;
      while(readVars==true){
        config>>dump;
        if(dump=="--EndAddVars--"){ readVars=false; continue;}
        OtherVars.push_back(new BDTVar);
        OtherVars.back()->name=dump;
      }
    }
  }
  config.close();
  
  std::cout<<"Config File read"<<std::endl;
  std::cout<<"--Initial Variables--"<<std::endl;
  for(int l=0;l<InitialVars.size();l++){std::cout<<InitialVars.at(l)->name<<std::endl;}
  std::cout<<"--AdditionalVariables--"<<std::endl;
  for(int l=0;l<OtherVars.size();l++){std::cout<<OtherVars.at(l)->name<<std::endl;}
  std::cout<<"--Method--"<<std::endl;
  std::cout<<MethodType<<std::endl;
  std::cout<<MethodString<<std::endl;
//done reading config

  
  Double_t KS=1.0;
  Double_t FOM=999.9;
  Double_t BestKS=0.0;
  Double_t BestFOM=0.0;
  Double_t SecBestKS=0.0;
  Double_t SecBestFOM=0.0;
  Double_t bufferKS=1.0;
  Double_t bufferFOM=999.9;

  std::vector<BDTVar*> BestVars;
  std::vector<BDTVar*> SecBestVars;
  std::vector<BDTVar*> WorstVars;
  std::vector<BDTVar*> UsedVars;
  std::vector<BDTVar*> UnusedVars;
  std::vector<BDTVar*> BestUnusedVars;

  //do initial training
  KS=1.0;
  FOM=999.9;
  for(int nn=0;nn<RepeatTrainingNTimes+1;nn++)
  {
    DoTraining(InitialVars,OtherVars,FOMType,FactoryString,PrepString,BWeight,CWeight,OWeight,BTreeName,CTreeName,OTreeName,MethodType,MethodString, particleNumber, &bufferFOM, &bufferKS,UseFixedTrainTestSplitting);
    //std::vector<BDTVar*> UsedVars,std::vector<BDTVar*> UnUsedVars,TString FOMType, TString FactoryString, TString PrepString, TString BWeight, TString CWeight, TString OWeight, TString BTreeName, TString CTreeName,TString OTreeName ,TString MethodType, TString MethodString, int particleNumber, Double_t* testFOM, Double_t* testKS,Int_t UseFixedTrainTestSplitting, TString PlotName="NONE" 

    if(bufferKS  < KS ){ KS  = bufferKS;  }
    if(bufferFOM < FOM){ FOM = bufferFOM; }
  }

  std::cout<<"Initial "<<RepeatTrainingNTimes+1<<" trainings KS, FOM "<<KS<<" "<<FOM<<std::endl;

  //check different Variable Combinations
  int nUsedVars=InitialVars.size();

  for(std::vector<BDTVar*>::iterator itVar=InitialVars.begin(); itVar!=InitialVars.end();++itVar){
    BDTVar* iVar=*itVar;
    BestVars.push_back(iVar);
    UsedVars.push_back(iVar);
  }

  for(std::vector<BDTVar*>::iterator itVar=OtherVars.begin(); itVar!=OtherVars.end();++itVar){
    BDTVar* iVar=*itVar;
    UnusedVars.push_back(iVar);
    WorstVars.push_back(iVar);
    BestUnusedVars.push_back(iVar);
  }

  if((KS < 1.0) and (KS > KSThreshold) and (FOM >= BestFOM))
  {
    BestFOM=FOM;
    BestKS=KS;
    BestVars.clear();
    for(std::vector<BDTVar*>::iterator itVar=InitialVars.begin(); itVar!=InitialVars.end();++itVar){
      BDTVar* iVar=*itVar;
      BestVars.push_back(iVar);
    }
  }

  std::cout<<"Initial KS , FOM "<<BestKS<<", "<<BestFOM<<std::endl<<std::endl;

  if(bool(FindBestVariables))
  {
    std::cout<<"Removing each Variable"<<std::endl;
    // remove worst Variable
    BDTVar* worstVar = new BDTVar;
    for(int k=0;k<nUsedVars;k++){
      BDTVar* testVar = new BDTVar;
      testVar->name=UsedVars.at(k)->name;
      UsedVars.erase(UsedVars.begin()+k);
      UnusedVars.push_back(testVar);
      
      KS=1.0;
      FOM=999.9;
      for(int nn=0;nn<RepeatTrainingNTimes+1;nn++){
      std::cout<<"Training Nr. "<<nn<<std::endl;  
      DoTraining(UsedVars,UnusedVars,FOMType,FactoryString,PrepString,BWeight,CWeight,OWeight,BTreeName,CTreeName,OTreeName,MethodType,MethodString, particleNumber, &bufferFOM, &bufferKS,UseFixedTrainTestSplitting);
      if(bufferKS<KS)KS=bufferKS;
      if(bufferFOM<FOM)FOM=bufferFOM;
      }
      std::cout<<"after "<<RepeatTrainingNTimes+1<<" trainigns KS, FOM "<<KS<<" "<<FOM<<std::endl;
      
      if(KS>KSThreshold and FOM>=BestFOM){
        BestFOM=FOM;
        BestKS=KS;
        BestVars.clear();
        for(int l=0;l<UsedVars.size();l++){
          BestVars.push_back(new BDTVar);
          BestVars.at(l)->name=UsedVars.at(l)->name; 
        }
        BestUnusedVars.clear();
        for(int l=0;l<UnusedVars.size();l++){
          BestUnusedVars.push_back(new BDTVar);
          BestUnusedVars.at(l)->name=UnusedVars.at(l)->name; 
        }
        
      }
//       std::cout<<FOM<<std::endl;
      if(KS>0.0 and FOM>SecBestFOM){
        SecBestFOM=FOM;
        SecBestKS=KS;
        SecBestVars.clear();
        for(int l=0;l<UsedVars.size();l++){
          SecBestVars.push_back(new BDTVar);
          SecBestVars.at(l)->name=UsedVars.at(l)->name; 
        }
        worstVar->name=testVar->name;
       }
      
       UsedVars.insert(UsedVars.begin()+k, testVar);
       UnusedVars.pop_back();
    }
    std::cout<<"worst Variable "<<worstVar->name<<std::endl;
    std::cout<<BestKS<<" "<<BestFOM<<std::endl;
    if(BestFOM==0.0){
      std::cout<<"here"<<std::endl;
      BestVars.clear();
      for(int l=0;l<SecBestVars.size();l++){
        BestVars.push_back(new BDTVar);
        BestVars.at(l)->name=SecBestVars.at(l)->name; 
      }
      BestUnusedVars.clear();
      for(int l=0;l<OtherVars.size();l++){
//         if(!OtherVars.at(l)->name.EqualTo(worstVar->name)){
          BestUnusedVars.push_back(new BDTVar);
          BestUnusedVars.back()->name=OtherVars.at(l)->name;
//         }
      }
      BestUnusedVars.push_back(new BDTVar);
      BestUnusedVars.back()->name=worstVar->name;
    }

    //add one more Variable to improve result
    BDTVar* addedVar=new BDTVar;
    BDTVar* testVar=new BDTVar;
    
    if(BestVars.size()<MaxVariablesInCombination || BestVars.size()<MinVariablesInCombination){
    
    std::cout<<"adding one variable"<<std::endl;
    UsedVars.clear();
    for(int k=0;k<SecBestVars.size();k++){
      UsedVars.push_back(new BDTVar);
      UsedVars.at(k)->name=SecBestVars.at(k)->name; 
    }
    OtherVars.push_back(new BDTVar);
    OtherVars.back()->name=worstVar->name;
    for(int k=0; k<OtherVars.size();k++){
      testVar->name=OtherVars.at(k)->name;
      UsedVars.push_back(testVar);
      OtherVars.erase(OtherVars.begin()+k);
      UnusedVars.clear();
      for(int l=0;l<OtherVars.size();l++){
        UnusedVars.push_back(new BDTVar);
        UnusedVars.at(l)->name=OtherVars.at(l)->name; 
      }
       
      KS=1.0;
      FOM=999.9;
      for(int nn=0;nn<RepeatTrainingNTimes+1;nn++){
      std::cout<<"Training Nr. "<<nn<<std::endl; 
      DoTraining(UsedVars,UnusedVars,FOMType,FactoryString,PrepString,BWeight,CWeight,OWeight,BTreeName,CTreeName,OTreeName,MethodType,MethodString, particleNumber, &bufferFOM, &bufferKS,UseFixedTrainTestSplitting);
      if(bufferKS<KS)KS=bufferKS;
      if(bufferFOM<FOM)FOM=bufferFOM;
      }
      std::cout<<"after "<<RepeatTrainingNTimes+1<<" trainigns KS, FOM "<<KS<<" "<<FOM<<std::endl;
      
      if(KS>KSThreshold and FOM>=ImprovementThreshold*BestFOM){
        BestFOM=FOM;
        BestKS=KS;
        BestVars.clear();
        for(int l=0;l<UsedVars.size();l++){
          BestVars.push_back(new BDTVar);
          BestVars.at(l)->name=UsedVars.at(l)->name; 
        }
        BestUnusedVars.clear();
        for(int l=0;l<UnusedVars.size();l++){
          BestUnusedVars.push_back(new BDTVar);
          BestUnusedVars.at(l)->name=UnusedVars.at(l)->name; 
        }
        addedVar->name=testVar->name;
        std::cout<<addedVar->name<<" "<<testVar->name<<std::endl;
      }
      
      OtherVars.insert(OtherVars.begin()+k, new BDTVar);
      OtherVars.at(k)->name=testVar->name;
      UsedVars.pop_back();
      
    }
//     std::cout<<addedVar->name<<" "<<testVar->name<<std::endl;

    std::cout<<"Variable added "<<addedVar->name<<std::endl;
    std::cout<<BestKS<<" "<<BestFOM<<std::endl;
    }
    //add another Variable to improve result
    
    
    //check if the previous step added a avriable
    bool dosecondTry=true;
    if(addedVar->name==""){
      std::cout<<"no improvement by adding variables"<<std::endl;
      dosecondTry=false;
    }
    if(BestVars.size()>=MaxVariablesInCombination){
      std::cout<<"maximal number of Variables to use reached "<<BestVars.size()<<std::endl;
      dosecondTry=false;
    }
    if(BestVars.size()<MinVariablesInCombination){
      std::cout<<" below the minimal number of Variables to use "<<BestVars.size()<<std::endl;
      dosecondTry=true;
    }
    if(dosecondTry){
    std::cout<<"add another"<<std::endl;
    BDTVar* secaddedVar=new BDTVar;
    UsedVars.clear();
    std::cout<<BestVars.size()<<std::endl;
    std::cout<<OtherVars.size()<<std::endl;
    
    

    for(int k=0;k<BestVars.size();k++){
      UsedVars.push_back(new BDTVar);
      UsedVars.at(k)->name=BestVars.at(k)->name; 
    }
    UnusedVars.clear();
//     std::cout<<"step 0"<<std::endl;

    for(int k=0;k<OtherVars.size();k++){
      std::cout<<"OtherVars "<<OtherVars.at(k)->name<<std::endl;

      if(!OtherVars.at(k)->name.EqualTo(addedVar->name)){
      UnusedVars.push_back(new BDTVar);
      UnusedVars.back()->name=OtherVars.at(k)->name;
      }
//       if(OtherVars.at(k)->name.EqualTo(addedVar->name))UnusedVars.pop_back(); 
    }
//     std::cout<<"step 1"<<std::endl;
    for(int k=0; k<UnusedVars.size();k++){
//       BDTVar* testVar;
      testVar->name=UnusedVars.at(k)->name;
      UsedVars.push_back(testVar);
      UnusedVars.erase(UnusedVars.begin()+k);
//     std::cout<<"step 2"<<std::endl;
         
      KS=1.0;
      FOM=999.9;
      for(int nn=0;nn<RepeatTrainingNTimes+1;nn++){
      std::cout<<"Training Nr. "<<nn<<std::endl;  
      DoTraining(UsedVars,UnusedVars,FOMType,FactoryString,PrepString,BWeight,CWeight,OWeight,BTreeName,CTreeName,OTreeName,MethodType,MethodString, particleNumber, &bufferFOM, &bufferKS,UseFixedTrainTestSplitting);
      if(bufferKS<KS)KS=bufferKS;
      if(bufferFOM<FOM)FOM=bufferFOM;
      }
      std::cout<<"after "<<RepeatTrainingNTimes+1<<" trainigns KS, FOM "<<KS<<" "<<FOM<<std::endl;

      
      if((KS>KSThreshold and FOM>=ImprovementThreshold*BestFOM) or (FOM>BestFOM and BestVars.size()<MinVariablesInCombination)){
        BestFOM=FOM;
        BestKS=KS;
        BestVars.clear();
        for(int l=0;l<UsedVars.size();l++){
          BestVars.push_back(new BDTVar);
          BestVars.at(l)->name=UsedVars.at(l)->name; 
        }
        WorstVars.clear();
        for(int l=0;l<UnusedVars.size();l++){
          WorstVars.push_back(new BDTVar);
          WorstVars.at(l)->name=UnusedVars.at(l)->name; 
        }
        BestUnusedVars.clear();
        for(int l=0;l<UnusedVars.size();l++){
          BestUnusedVars.push_back(new BDTVar);
          BestUnusedVars.at(l)->name=UnusedVars.at(l)->name; 
        }
        secaddedVar->name=testVar->name;
        std::cout<<secaddedVar->name<<" "<<testVar->name<<std::endl;

      }
      
      UnusedVars.insert(UnusedVars.begin()+k, new BDTVar);
      UnusedVars.at(k)->name=testVar->name;
      UsedVars.pop_back();
      
    }
    std::cout<<"second Variable added "<<secaddedVar->name<<std::endl;
    std::cout<<BestKS<<" "<<BestFOM<<std::endl;
    }
    else{
      std::cout<<"return to original Vars"<<std::endl;
//       BestVars.clear();
//       for(int l=0;l<SecBestVars.size();l++){
//         BestVars.push_back(new BDTVar);
//         BestVars.at(l)->name=SecBestVars.at(l)->name; 
//       }
      WorstVars.clear();
      for(int l=0;l<OtherVars.size();l++){
        WorstVars.push_back(new BDTVar);
        WorstVars.at(l)->name=OtherVars.at(l)->name; 
      }  
    }

    std::cout<<"Best FOM "<<BestFOM<<" "<<BestKS<<std::endl;

    for(int k=0;k<BestVars.size();k++){
      std::cout<<BestVars.at(k)->name<<std::endl; 
    }

//    TString ROCFileName="ROC_Particle";
//    ROCFileName+=particleNumber;
//    ROCFileName+="_Iteration";
//    ROCFileName+=Iteration;

    KS=1.0;
    FOM=999.9;
    for(int nn=0;nn<RepeatTrainingNTimes+1;nn++)
    {
      DoTraining(BestVars,WorstVars,FOMType,FactoryString,PrepString,BWeight,CWeight,OWeight,BTreeName,CTreeName,OTreeName,MethodType,MethodString, particleNumber, &bufferFOM, &bufferKS,UseFixedTrainTestSplitting);
      if(bufferKS<KS)KS=bufferKS;
      if(bufferFOM<FOM)FOM=bufferFOM;
    }

    std::cout<<"Final "<<RepeatTrainingNTimes+1<<" trainigns KS, FOM "<<KS<<" "<<FOM<<std::endl;

    BestFOM = FOM;
    BestKS = KS;
  }

  std::ofstream result("ParticleResult.txt", std::ofstream::trunc);
  result<<"BestFOM "<<BestFOM<<std::endl;
  result<<"KSScore "<<BestKS<<std::endl;
  result<<"MethodString "<<MethodString<<std::endl;
//   result<<BestVars.size()<<std::endl;
//   result<<BestUnusedVars.size()<<std::endl;
  for(int k=0;k<BestVars.size();k++){
    std::cout<<"UsedVar "<<BestVars.at(k)->name<<std::endl; 
    result<<"UsedVar "<<BestVars.at(k)->name<<std::endl; 
  }
//   result<<"UnusedVars"<<std::endl;
  for(int k=0;k<BestUnusedVars.size();k++){
    std::cout<<"UnusedVars "<<BestUnusedVars.at(k)->name<<std::endl;
    result<<"UnusedVars "<<BestUnusedVars.at(k)->name<<std::endl; 
  }
  result.close();

  TimeTotal = timerTotal->RealTime();

  delete timerTotal;

  std::cout<<"*******************************************************"<<std::endl;
  std::cout<<"* Total Time: "<<TimeTotal<<std::endl;
  std::cout<<"* Tree processing time : "<<TimeTreePreparation<< " "<<TimeTreePreparation/TimeTotal*100.0<<"%"<<std::endl;
  std::cout<<"* Training time : "<<TimeTraining<< " "<<TimeTraining/TimeTotal*100.0<<"%"<<std::endl;
  std::cout<<"* Testing time : "<<TimeTesting<< " "<<TimeTesting/TimeTotal*100.0<<"%"<<std::endl;
  std::cout<<"*******************************************************"<<std::endl;
}


# ifndef __CINT__
int main()
{
  Particle();
  return 0;
}
# endif
