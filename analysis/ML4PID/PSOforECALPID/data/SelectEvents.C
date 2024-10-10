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

void SelectEvents(TString filename, TString detector, int NEvents) {
  std::cout<<"Loading file : "<<filename<<std::endl;
  TFile *f = new TFile(filename+".root","update");
  TTree *T = (TTree*)f->Get("ntp"); 
  
  Float_t ecal_interaction = 0;
  T->SetBranchAddress("ecal_interaction",&ecal_interaction);
  Float_t hcal_interaction = 0;
  T->SetBranchAddress("hcal_interaction",&hcal_interaction);
  Float_t total_interaction = 0;
  T->SetBranchAddress("total_interaction",&total_interaction);
  
  Long64_t nentries = T->GetEntries();
  
  TFile newfile(filename+"_"+detector+".root", "recreate");
  auto newtree = T->CloneTree(0);
  
  int eventcount = 0;

  for (int i = 0; i<nentries/3-3; i++) {
    T->GetEntry(3*i);           
    //cout<<"(e, h, t) = ("<<ecal_interaction<<", "<<hcal_interaction<<", "<<total_interaction<<")"<<endl;
    if((detector=="ecal") and (ecal_interaction==1)) {
      newtree->Fill();
      eventcount += 1;
    }
    if((detector=="hcal") and(hcal_interaction==1)) {
      newtree->Fill();
      eventcount += 1;
    }
    if((detector=="total") and(total_interaction==1)) {
      newtree->Fill();
      eventcount += 1;
    }
   
    if(eventcount == NEvents) {
      cout<<"Reached "<<NEvents<<" events."<<endl;
      break;
    }    
    //e_inter->Clear();
    //h_inter->Clear();
    //t_inter->Clear();
  }

  for (int i = 0; i<nentries/3-3; i++) {
    T->GetEntry(3*i+1);
    //cout<<"(e, h, t) = ("<<ecal_interaction<<", "<<hcal_interaction<<", "<<total_interaction<<")"<<endl;
    if(eventcount == NEvents) {
      break;
    }
    if((detector=="ecal") and (ecal_interaction==1)) {
      newtree->Fill();
      eventcount += 1;
    }
    if((detector=="hcal") and(hcal_interaction==1)) {
      newtree->Fill();
      eventcount += 1;
    }
    if((detector=="total") and(total_interaction==1)) {
      newtree->Fill();
      eventcount += 1;
    }

    if(eventcount == NEvents) {
      cout<<"Reached "<<NEvents<<" events."<<endl;
      break;
    }
    //e_inter->Clear();
    //h_inter->Clear();
    //t_inter->Clear();
  }

  for (int i = 0; i<nentries/3-3; i++) {
    T->GetEntry(3*i+2);
    //cout<<"(e, h, t) = ("<<ecal_interaction<<", "<<hcal_interaction<<", "<<total_interaction<<")"<<endl;               
    if(eventcount == NEvents) {
      break;
    }
    if((detector=="ecal") and (ecal_interaction==1)) {
      newtree->Fill();
      eventcount += 1;
    }
    if((detector=="hcal") and(hcal_interaction==1)) {
      newtree->Fill();
      eventcount += 1;
    }
    if((detector=="total") and(total_interaction==1)) {
      newtree->Fill();
      eventcount += 1;
    }

    if(eventcount == NEvents) {
      cout<<"Reached "<<NEvents<<" events."<<endl;
      break;
    }
    //e_inter->Clear();
    //h_inter->Clear();
    //t_inter->Clear();
  }
  //T->Print(); 
  newtree->Write(); 
  delete f;

  if(eventcount != NEvents) cout<<"Didn't get to "<<NEvents<<" events :( "<<endl;
  
}
