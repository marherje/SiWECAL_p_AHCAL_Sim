#include <TPaveStats.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <TFitResult.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystemFile.h"
#include "TGraph.h"
#include "TAxis.h"

void AddWeightBranch(TString filename) {

  std::cout<<"Loading file "<<std::endl;
  TFile *f = new TFile(filename,"update");
  TTree *T = (TTree*)f->Get("ntp"); 
  //float px,py; 
  //float pt;
  //TBranch *bpt = T->Branch("pt",&pt,"pt/F");
  float peso;
  TBranch *Weight = T->Branch("Weight",&peso,"Weight/F");
  //T->SetBranchAddress("px",&px);
  //T->SetBranchAddress("py",&py);
  Long64_t nentries = T->GetEntries();
  for (Long64_t i=0;i<nentries;i++) {
    T->GetEntry(i);
    //pt = TMath::Sqrt(px*px+py*py);
    //bpt->Fill();
    peso=1.0;
    Weight->Fill();
} 
  std::cout<<"Weight succesfully added"<<std::endl;
  T->Print(); 
  T->Write(); 
  delete f;

}
