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
#include "TH1.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystemFile.h"
#include "style/Style.C"
#include "style/Labels.C"
#include "TVirtualFitter.h"

float ptosigma(float p){
  float sigma;
  if(p==0)sigma=5;
  else if(p<0.000007)sigma=5;
  else if(p<0.00007)sigma=4.5;
  else if(p<0.00011)sigma=4;
  else if(p<0.0005)sigma=3.5;
  else if(p<0.0027)sigma=3;
  else if(p<0.0125)sigma=2.5;
  else if(p<0.05)sigma=2;
  else if(p<0.134)sigma=1.5;
  else if(p<0.317)sigma=1;
  else if(p<0.617)sigma=0.5;
  else sigma=0.001;

  return sigma;
}

void lcfiplusdEdxKS(const int samples,TString cat){

  // Data:                                                                                                                                                                    
  std::cout<<"Init Reading"<<std::endl;

  TString foldername="/lustre/ific.uv.es/prj/ific/flc/ntuples-2022/250/dEdx_Samples/";
  TString quark_string[3]={"b_quark","c_quark","light_quark"};
  //std::vector<string> quark_vector;
  //std::vector<string> sample_vector;
  //TString sample_string[size]={};
  /*
  if(samples==2){
    sample_string.push_back("00015161_eL_pR");
    sample_string.push_back("00015164_eR_pL");
  }
  */
  
  TString sample_string[5]={};
  if(samples==2){
    sample_string[0]="00015161_eL_pR";
    sample_string[1]="00015164_eR_pL";
  }
  else if(samples==3){
    sample_string[0]="00015162_eL_pR";
    sample_string[1]="00015165_eR_pL";
    sample_string[2]="00015275_eR_pL";
  }
  else if(samples==5){
    sample_string[0]="00015161_eL_pR";
    sample_string[1]="00015162_eL_pR";
    sample_string[2]="00015164_eR_pL";
    sample_string[3]="00015165_eR_pL";
    sample_string[4]="00015275_eR_pL";
  }

  //TString sample_string[5]={"merged","merged","merged","merged","merged"};
  TString filename[3][5];
  //std::vector<string> filename_b;
  //std::vector<string> filename_c;
  //std::vector<string> filename_uds;

  //merged_${quark}_${nameaux}_eL_pR.root
  TFile *f[3][5];
  std::cout<<"Files loaded: "<<std::endl;

  
  
  // Categories:
  // cat A = nvtx==0;
  // cat B = nvtx==1&&nvtxall==1;
  // cat C = nvtx==1&&nvtxall==2;
  // cat D = nvtx>=2;

  // Things we need for reading data
  // We read once, then we store in histos.
  TTree *tree;
  // Branches
  TBranch *nvtx;
  TBranch *nvtxall; //control
  // Good old variable
  TBranch *d0bprob2;  
  // New variables
  TBranch *dEdxRatioPionOverKaonPri;
  TBranch *dEdxRatioKaonOverProtonPri;
  TBranch *dEdxRatioPionOverProtonPri;
  TBranch *dEdxRatioPionOverKaonSec;
  TBranch *dEdxRatioKaonOverProtonSec;
  TBranch *dEdxRatioPionOverProtonSec;
  TBranch *dEdxNKaonPri;
  TBranch *dEdxNPionPri;
  TBranch *dEdxNProtonPri;
  TBranch *dEdxNKaonSec;
  TBranch *dEdxNPionSec;
  TBranch *dEdxNProtonSec;

  // Addresses                                                                                                     
  Float_t nvtx_a;
  Float_t nvtxall_a; //control
  Float_t d0bprob2_a;
  Float_t dEdxRatioPionOverKaonPri_a;
  Float_t dEdxRatioKaonOverProtonPri_a;
  Float_t dEdxRatioPionOverProtonPri_a;
  Float_t dEdxRatioPionOverKaonSec_a;
  Float_t dEdxRatioKaonOverProtonSec_a;
  Float_t dEdxRatioPionOverProtonSec_a;
  Float_t dEdxNKaonSec_a;
  Float_t dEdxNPionSec_a;
  Float_t dEdxNProtonSec_a;
  Float_t dEdxNKaonPri_a;
  Float_t dEdxNPionPri_a;
  Float_t dEdxNProtonPri_a;


  TH1F *d0bprob2_histo[3][5]; //[flavour][Samples]
  //TCanvas *d0bprob2_canvas = new TCanvas("d0bprob2 canvas","d0bprob2 canvas",1000,1000);
  TH1F *dEdxRatioPionOverKaonPri_histo[3][5];
  TH1F *dEdxRatioKaonOverProtonPri_histo[3][5];
  TH1F *dEdxRatioPionOverProtonPri_histo[3][5];
  //TCanvas *dEdxRatioPionOverKaonPri_canvas = new TCanvas("dEdxRatioPionOverKaonPri canvas","dEdxRatioPionOverKaonPri canvas",1000,1000);
  TH1F *dEdxRatioPionOverKaonSec_histo[3][5];
  TH1F *dEdxRatioKaonOverProtonSec_histo[3][5];
  TH1F *dEdxRatioPionOverProtonSec_histo[3][5];
  //TCanvas *dEdxRatioPionOverKaonSec_canvas = new TCanvas("dEdxRatioPionOverKaonSec canvas","dEdxRatioPionOverKaonSec canvas",1000,1000);
  TH1F *dEdxNKaonPri_histo[3][5];
  TH1F *dEdxNPionPri_histo[3][5];
  TH1F *dEdxNProtonPri_histo[3][5];
  TH1F *dEdxNKaonSec_histo[3][5];
  TH1F *dEdxNPionSec_histo[3][5];
  TH1F *dEdxNProtonSec_histo[3][5];

  for(int i=0;i<3;i++){
    for(int j=0;j<samples;j++){
      //merged_${quark}_${nameaux}_eL_pR.root
      filename[i][j]=foldername+"merged_"+quark_string[i]+"_"+sample_string[j]+".root";
      
      //filename[i][j]="../../data/MakeNTuples250_dEdx/"+quark_string[i]+"_"+sample_string[j]+".root";

      //eL_pR_b_quark_merged.root
      f[i][j]=new TFile(filename[i][j]);
      std::cout<<filename[i][j]<<std::endl;
            
      // Load tree
      tree = f[i][j]->Get<TTree>("ntp");
      // Load branches
      nvtx = tree->GetBranch("nvtx");
      nvtxall = tree->GetBranch("nvtxall");
      d0bprob2 = tree->GetBranch("d0bprob2");
      dEdxRatioPionOverKaonPri = tree->GetBranch("dEdxRatioPionOverKaonPri");
      dEdxRatioKaonOverProtonPri = tree->GetBranch("dEdxRatioKaonOverProtonPri");
      dEdxRatioPionOverProtonPri = tree->GetBranch("dEdxRatioPionOverProtonPri");
      dEdxRatioPionOverKaonSec = tree->GetBranch("dEdxRatioPionOverKaonSec");
      dEdxRatioKaonOverProtonSec = tree->GetBranch("dEdxRatioKaonOverProtonSec");
      dEdxRatioPionOverProtonSec = tree->GetBranch("dEdxRatioPionOverProtonSec");
      dEdxNKaonPri = tree->GetBranch("dEdxNKaonPri");
      dEdxNPionPri = tree->GetBranch("dEdxNPionPri");
      dEdxNProtonPri = tree->GetBranch("dEdxNProtonPri");
      dEdxNKaonSec = tree->GetBranch("dEdxNKaonSec");
      dEdxNPionSec = tree->GetBranch("dEdxNPionSec");
      dEdxNProtonSec = tree->GetBranch("dEdxNProtonSec");
      // Set Addresses                                                                                                                                             
      nvtx->SetAddress(&nvtx_a);
      nvtxall->SetAddress(&nvtxall_a);
      d0bprob2->SetAddress(&d0bprob2_a);
      dEdxRatioPionOverKaonPri->SetAddress(&dEdxRatioPionOverKaonPri_a);
      dEdxRatioKaonOverProtonPri->SetAddress(&dEdxRatioKaonOverProtonPri_a);
      dEdxRatioPionOverProtonPri->SetAddress(&dEdxRatioPionOverProtonPri_a);
      dEdxRatioPionOverKaonSec->SetAddress(&dEdxRatioPionOverKaonSec_a);
      dEdxRatioKaonOverProtonSec->SetAddress(&dEdxRatioKaonOverProtonSec_a);
      dEdxRatioPionOverProtonSec->SetAddress(&dEdxRatioPionOverProtonSec_a);
      dEdxNKaonPri->SetAddress(&dEdxNKaonPri_a);
      dEdxNPionPri->SetAddress(&dEdxNPionPri_a);
      dEdxNProtonPri->SetAddress(&dEdxNProtonPri_a);
      dEdxNKaonSec->SetAddress(&dEdxNKaonSec_a);
      dEdxNPionSec->SetAddress(&dEdxNPionSec_a);
      dEdxNProtonSec->SetAddress(&dEdxNProtonSec_a);
      
      int totalentries;
      totalentries = d0bprob2->GetEntries();
      std::cout<<"File entries: "<<totalentries<<std::endl;
      
      d0bprob2_histo[i][j] = new TH1F("d0bprob2 "+quark_string[i]+" "+sample_string[j],"d0bprob2 "+quark_string[i]+" "+sample_string[j],50,0,1);
      dEdxRatioPionOverKaonPri_histo[i][j] = new TH1F("dEdxRatioPionOverKaonPri "+quark_string[i]+" "+sample_string[j],"dEdxRatioPionOverKaonPri "+quark_string[i]+" "+sample_string[j],23,-2,20);
      dEdxRatioKaonOverProtonPri_histo[i][j] = new TH1F("dEdxRatioKaonOverProtonPri "+quark_string[i]+" "+sample_string[j],"dEdxRatioKaonOverProtonPri "+quark_string[i]+" "+sample_string[j],23,-2,20);
      dEdxRatioPionOverProtonPri_histo[i][j] = new TH1F("dEdxRatioPionOverProtonPri "+quark_string[i]+" "+sample_string[j],"dEdxRatioPionOverProtonPri "+quark_string[i]+" "+sample_string[j],23,-2,20);
      dEdxRatioPionOverKaonSec_histo[i][j] = new TH1F("dEdxRatioPionOverKaonSec "+quark_string[i]+" "+sample_string[j],"dEdxRatioPionOverKaonSec "+quark_string[i]+" "+sample_string[j],23,-2,20);
      dEdxRatioKaonOverProtonSec_histo[i][j] = new TH1F("dEdxRatioKaonOverProtonSec "+quark_string[i]+" "+sample_string[j],"dEdxRatioKaonOverProtonSec "+quark_string[i]+" "+sample_string[j],23,-2,20);
      dEdxRatioPionOverProtonSec_histo[i][j] = new TH1F("dEdxRatioPionOverProtonSec "+quark_string[i]+" "+sample_string[j],"dEdxRatioPionOverProtonSec "+quark_string[i]+" "+sample_string[j],23,-2,20);
      dEdxNKaonPri_histo[i][j] = new TH1F("dEdxNKaonPri "+quark_string[i]+" "+sample_string[j],"dEdxNKaonPri "+quark_string[i]+" "+sample_string[j],20,0,20);
      dEdxNPionPri_histo[i][j] = new TH1F("dEdxNPionPri "+quark_string[i]+" "+sample_string[j],"dEdxNPionPri "+quark_string[i]+" "+sample_string[j],20,0,20);
      dEdxNProtonPri_histo[i][j] = new TH1F("dEdxNProtonPri "+quark_string[i]+" "+sample_string[j],"dEdxNProtonPri "+quark_string[i]+" "+sample_string[j],20,0,20);
      dEdxNKaonSec_histo[i][j] = new TH1F("dEdxNKaonSec "+quark_string[i]+" "+sample_string[j],"dEdxNKaonSec "+quark_string[i]+" "+sample_string[j],20,0,20);
      dEdxNPionSec_histo[i][j] = new TH1F("dEdxNPionSec "+quark_string[i]+" "+sample_string[j],"dEdxNPionSec "+quark_string[i]+" "+sample_string[j],20,0,20);
      dEdxNProtonSec_histo[i][j] = new TH1F("dEdxNProtonSec "+quark_string[i]+" "+sample_string[j],"dEdxNProtonSec "+quark_string[i]+" "+sample_string[j],20,0,20);

      for(int ientry=0;ientry<totalentries;ientry++){
	nvtx->GetEntry(ientry);
	nvtxall->GetEntry(ientry);
	
	d0bprob2->GetEntry(ientry);
	dEdxRatioPionOverKaonPri->GetEntry(ientry);
	dEdxRatioKaonOverProtonPri->GetEntry(ientry);
	dEdxRatioPionOverProtonPri->GetEntry(ientry);
	dEdxRatioPionOverKaonSec->GetEntry(ientry);
	dEdxRatioKaonOverProtonSec->GetEntry(ientry);
	dEdxRatioPionOverProtonSec->GetEntry(ientry);
	dEdxNKaonPri->GetEntry(ientry);
	dEdxNPionPri->GetEntry(ientry);
	dEdxNProtonPri->GetEntry(ientry);
	dEdxNKaonSec->GetEntry(ientry);
	dEdxNPionSec->GetEntry(ientry);
	dEdxNProtonSec->GetEntry(ientry);

	// Fill the corresponding category
	if(cat == "all"){
	  d0bprob2_histo[i][j]->Fill(d0bprob2_a);
	  dEdxRatioPionOverKaonPri_histo[i][j]->Fill(dEdxRatioPionOverKaonPri_a);
	  dEdxRatioKaonOverProtonPri_histo[i][j]->Fill(dEdxRatioKaonOverProtonPri_a);
	  dEdxRatioPionOverProtonPri_histo[i][j]->Fill(dEdxRatioPionOverProtonPri_a);
	  dEdxRatioPionOverKaonSec_histo[i][j]->Fill(dEdxRatioPionOverKaonSec_a);
	  dEdxRatioKaonOverProtonSec_histo[i][j]->Fill(dEdxRatioKaonOverProtonSec_a);
	  dEdxRatioPionOverProtonSec_histo[i][j]->Fill(dEdxRatioPionOverProtonSec_a);
	  dEdxNKaonPri_histo[i][j]->Fill(dEdxNKaonPri_a);
	  dEdxNPionPri_histo[i][j]->Fill(dEdxNPionPri_a);
	  dEdxNProtonPri_histo[i][j]->Fill(dEdxNProtonPri_a);
	  dEdxNKaonSec_histo[i][j]->Fill(dEdxNKaonSec_a);
	  dEdxNPionSec_histo[i][j]->Fill(dEdxNPionSec_a);
	  dEdxNProtonSec_histo[i][j]->Fill(dEdxNProtonSec_a);
	}
	else if(cat == "catA"){
	  if(nvtx_a==0){
	    d0bprob2_histo[i][j]->Fill(d0bprob2_a);
	    dEdxRatioPionOverKaonPri_histo[i][j]->Fill(dEdxRatioPionOverKaonPri_a);
	    dEdxRatioKaonOverProtonPri_histo[i][j]->Fill(dEdxRatioKaonOverProtonPri_a);
	    dEdxRatioPionOverProtonPri_histo[i][j]->Fill(dEdxRatioPionOverProtonPri_a);
	    dEdxRatioPionOverKaonSec_histo[i][j]->Fill(dEdxRatioPionOverKaonSec_a);
	    dEdxRatioKaonOverProtonSec_histo[i][j]->Fill(dEdxRatioKaonOverProtonSec_a);
	    dEdxRatioPionOverProtonSec_histo[i][j]->Fill(dEdxRatioPionOverProtonSec_a);
	    dEdxNKaonPri_histo[i][j]->Fill(dEdxNKaonPri_a);
	    dEdxNPionPri_histo[i][j]->Fill(dEdxNPionPri_a);
	    dEdxNProtonPri_histo[i][j]->Fill(dEdxNProtonPri_a);
	    dEdxNKaonSec_histo[i][j]->Fill(dEdxNKaonSec_a);
	    dEdxNPionSec_histo[i][j]->Fill(dEdxNPionSec_a);
	    dEdxNProtonSec_histo[i][j]->Fill(dEdxNProtonSec_a);
	  }
	}
	else if(cat == "catB"){
	  if(nvtx_a==1&&nvtxall_a==1){
	    d0bprob2_histo[i][j]->Fill(d0bprob2_a);
            dEdxRatioPionOverKaonPri_histo[i][j]->Fill(dEdxRatioPionOverKaonPri_a);
	    dEdxRatioKaonOverProtonPri_histo[i][j]->Fill(dEdxRatioKaonOverProtonPri_a);
	    dEdxRatioPionOverProtonPri_histo[i][j]->Fill(dEdxRatioPionOverProtonPri_a);
	    dEdxRatioPionOverKaonSec_histo[i][j]->Fill(dEdxRatioPionOverKaonSec_a);
	    dEdxRatioKaonOverProtonSec_histo[i][j]->Fill(dEdxRatioKaonOverProtonSec_a);
	    dEdxRatioPionOverProtonSec_histo[i][j]->Fill(dEdxRatioPionOverProtonSec_a);
	    dEdxNKaonPri_histo[i][j]->Fill(dEdxNKaonPri_a);
	    dEdxNPionPri_histo[i][j]->Fill(dEdxNPionPri_a);
	    dEdxNProtonPri_histo[i][j]->Fill(dEdxNProtonPri_a);
	    dEdxNKaonSec_histo[i][j]->Fill(dEdxNKaonSec_a);
	    dEdxNPionSec_histo[i][j]->Fill(dEdxNPionSec_a);
	    dEdxNProtonSec_histo[i][j]->Fill(dEdxNProtonSec_a);
	  }
	}
	else if(cat == "catC"){
	  if(nvtx_a==1&&nvtxall_a==2){
	    d0bprob2_histo[i][j]->Fill(d0bprob2_a);
	    dEdxRatioPionOverKaonPri_histo[i][j]->Fill(dEdxRatioPionOverKaonPri_a);
	    dEdxRatioKaonOverProtonPri_histo[i][j]->Fill(dEdxRatioKaonOverProtonPri_a);
	    dEdxRatioPionOverProtonPri_histo[i][j]->Fill(dEdxRatioPionOverProtonPri_a);
	    dEdxRatioPionOverKaonSec_histo[i][j]->Fill(dEdxRatioPionOverKaonSec_a);
	    dEdxRatioKaonOverProtonSec_histo[i][j]->Fill(dEdxRatioKaonOverProtonSec_a);
	    dEdxRatioPionOverProtonSec_histo[i][j]->Fill(dEdxRatioPionOverProtonSec_a);
	    dEdxNKaonPri_histo[i][j]->Fill(dEdxNKaonPri_a);
	    dEdxNPionPri_histo[i][j]->Fill(dEdxNPionPri_a);
	    dEdxNProtonPri_histo[i][j]->Fill(dEdxNProtonPri_a);
	    dEdxNKaonSec_histo[i][j]->Fill(dEdxNKaonSec_a);
	    dEdxNPionSec_histo[i][j]->Fill(dEdxNPionSec_a);
	    dEdxNProtonSec_histo[i][j]->Fill(dEdxNProtonSec_a);
	  }
	}
	else if(cat == "catD"){
	  if(nvtx_a==2){
	    d0bprob2_histo[i][j]->Fill(d0bprob2_a);
	    dEdxRatioPionOverKaonPri_histo[i][j]->Fill(dEdxRatioPionOverKaonPri_a);
	    dEdxRatioKaonOverProtonPri_histo[i][j]->Fill(dEdxRatioKaonOverProtonPri_a);
	    dEdxRatioPionOverProtonPri_histo[i][j]->Fill(dEdxRatioPionOverProtonPri_a);
	    dEdxRatioPionOverKaonSec_histo[i][j]->Fill(dEdxRatioPionOverKaonSec_a);
	    dEdxRatioKaonOverProtonSec_histo[i][j]->Fill(dEdxRatioKaonOverProtonSec_a);
	    dEdxRatioPionOverProtonSec_histo[i][j]->Fill(dEdxRatioPionOverProtonSec_a);
	    dEdxNKaonPri_histo[i][j]->Fill(dEdxNKaonPri_a);
	    dEdxNPionPri_histo[i][j]->Fill(dEdxNPionPri_a);
	    dEdxNProtonPri_histo[i][j]->Fill(dEdxNProtonPri_a);
	    dEdxNKaonSec_histo[i][j]->Fill(dEdxNKaonSec_a);
	    dEdxNPionSec_histo[i][j]->Fill(dEdxNPionSec_a);
	    dEdxNProtonSec_histo[i][j]->Fill(dEdxNProtonSec_a);
	  }
	}
      }//entries
	
      int entriesincat = d0bprob2_histo[i][j]->GetEntries();
      std::cout<<"Entries in d0bprob2 in cat. "<<cat<<": "<<entriesincat<<std::endl;
      Double_t factor=1.;
      d0bprob2_histo[i][j]->Scale(factor/entriesincat);
      dEdxRatioPionOverKaonPri_histo[i][j]->Scale(factor/entriesincat);
      dEdxRatioKaonOverProtonPri_histo[i][j]->Scale(factor/entriesincat);
      dEdxRatioPionOverProtonPri_histo[i][j]->Scale(factor/entriesincat);
      dEdxRatioPionOverKaonSec_histo[i][j]->Scale(factor/entriesincat);
      dEdxRatioKaonOverProtonSec_histo[i][j]->Scale(factor/entriesincat);
      dEdxRatioPionOverProtonSec_histo[i][j]->Scale(factor/entriesincat);
      dEdxNKaonPri_histo[i][j]->Scale(factor/entriesincat);
      dEdxNPionPri_histo[i][j]->Scale(factor/entriesincat);
      dEdxNProtonPri_histo[i][j]->Scale(factor/entriesincat);
      dEdxNKaonSec_histo[i][j]->Scale(factor/entriesincat);
      dEdxNPionSec_histo[i][j]->Scale(factor/entriesincat);
      dEdxNProtonSec_histo[i][j]->Scale(factor/entriesincat);
    
    }//j
  }//i
  
  // KS Test
  // 3 jxj matrices (histos)
  float ks_d0bprob2[3][5][5]={};
  Float_t ks_dEdxRatioPionOverKaonPri[3][5][5]={};
  Float_t ks_dEdxRatioPionOverKaonSec[3][5][5]={};
  Float_t ks_dEdxNKaonPri[3][5][5]={};
  Float_t ks_dEdxNKaonSec[3][5][5]={};
  // 2D histo for the cross-comparison
  TH2F *KS_d0bprob2[3];
  //TH2F *AD_d0bprob2[3];
  TH2F *KS_dEdxRatioPionOverKaonPri[3];
  TH2F *KS_dEdxRatioPionOverKaonSec[3];
  TH2F *KS_dEdxNKaonPri[3];
  TH2F *KS_dEdxNKaonSec[3];
  
  for(int i=0;i<3;i++){
    for(int j=0;j<samples;j++){ //Samples
      for(int k=0;k<samples;k++){ //Samples
	//ks_d0bprob2[i][j][k]=d0bprob2_histo[i][j]->KolmogorovTest(d0bprob2_histo[i][k], "X");
	//ks_dEdxRatioPionOverKaonPri[i][j][k]=dEdxRatioPionOverKaonPri_histo[i][j]->KolmogorovTest(dEdxRatioPionOverKaonPri_histo[i][k], "X");
	//ks_dEdxRatioPionOverKaonSec[i][j][k]=dEdxRatioPionOverKaonSec_histo[i][j]->KolmogorovTest(dEdxRatioPionOverKaonSec_histo[i][k], "X");
	//ks_dEdxNKaonPri[i][j][k]=dEdxNKaonPri_histo[i][j]->KolmogorovTest(dEdxNKaonPri_histo[i][k], "X");
	//ks_dEdxNKaonSec[i][j][k]=dEdxNKaonSec_histo[i][j]->KolmogorovTest(dEdxNKaonSec_histo[i][k], "X");
	
	int av_limit=30;
	Float_t KS_d0bprob2_av=0;
	Float_t KS_d0bprob2_aux=0;
	Float_t KS_KoP_p_av=0;
	Float_t KS_KoP_p_aux=0;
	Float_t KS_KoP_s_av=0;
	Float_t KS_KoP_s_aux=0;
	Float_t KS_NK_p_av=0;
	Float_t KS_NK_p_aux=0;
	Float_t KS_NK_s_av=0;
	Float_t KS_NK_s_aux=0;
	for(int z=0;z<av_limit;z++){
	  KS_d0bprob2_aux=d0bprob2_histo[i][j]->KolmogorovTest(d0bprob2_histo[i][k], "X");
	  KS_d0bprob2_av+=KS_d0bprob2_aux;
	  KS_KoP_p_aux=dEdxRatioPionOverKaonPri_histo[i][j]->KolmogorovTest(dEdxRatioPionOverKaonPri_histo[i][k], "X");
          KS_KoP_p_av+=KS_KoP_p_aux;
	  KS_KoP_s_aux=dEdxRatioPionOverKaonSec_histo[i][j]->KolmogorovTest(dEdxRatioPionOverKaonSec_histo[i][k], "X");
          KS_KoP_s_av+=KS_KoP_s_aux;
	  KS_NK_p_aux=dEdxNKaonPri_histo[i][j]->KolmogorovTest(dEdxNKaonPri_histo[i][k], "X");
          KS_NK_p_av+=KS_NK_p_aux;
	  KS_NK_s_aux=dEdxNKaonSec_histo[i][j]->KolmogorovTest(dEdxNKaonSec_histo[i][k], "X");
          KS_NK_s_av+=KS_NK_s_aux;
	}
	ks_d0bprob2[i][j][k]=KS_d0bprob2_av/av_limit;
	ks_dEdxRatioPionOverKaonPri[i][j][k]=KS_KoP_p_av/av_limit;
	ks_dEdxRatioPionOverKaonSec[i][j][k]=KS_KoP_s_av/av_limit;
	ks_dEdxNKaonPri[i][j][k]=KS_NK_p_av/av_limit;
	ks_dEdxNKaonSec[i][j][k]=KS_NK_s_av/av_limit;
	//cout<<"KS d0bprob2: "<<ks_d0bprob2[i][j][k]<<endl;
	//cout<<"KS PoKpri: "<<ks_dEdxRatioPionOverKaonPri[i][j][k]<<endl;
	//cout<<"KS PoKsec: "<<ks_dEdxRatioPionOverKaonSec[i][j][k]<<endl;
	//cout<<"KS NKpri: "<<ks_dEdxNKaonPri[i][j][k]<<endl; 
	//cout<<"KS NKsec: "<<ks_dEdxNKaonSec[i][j][k]<<endl;
      }  
    }
  }  
  
  float hlimit=samples-0.5;
  for(int i=0;i<3;i++){
    KS_d0bprob2[i] = new TH2F("KS d0bprob2 "+quark_string[i]+" flavour", "KS d0bprob2 "+quark_string[i]+" flavour",samples,-0.5,hlimit,samples,-0.5,hlimit);
    KS_dEdxRatioPionOverKaonPri[i] = new TH2F("KS dEdxRatioPionOverKaonPri "+quark_string[i]+" flavour", "KS dEdxRatioPionOverKaonPri "+quark_string[i]+" flavour",samples,-0.5,hlimit,samples,-0.5,hlimit);
    KS_dEdxRatioPionOverKaonSec[i] = new TH2F("KS dEdxRatioPionOverKaonSec "+quark_string[i]+" flavour", "KS dEdxRatioPionOverKaonSec "+quark_string[i]+" flavour",samples,-0.5,hlimit,samples,-0.5,hlimit);
    KS_dEdxNKaonPri[i] = new TH2F("KS dEdxNKaonPri "+quark_string[i]+" flavour", "KS dEdxNKaonPri "+quark_string[i]+" flavour",samples,-0.5,hlimit,samples,-0.5,hlimit);
    KS_dEdxNKaonSec[i] = new TH2F("KS dEdxNKaonSec "+quark_string[i]+" flavour", "KS dEdxNKaonSec "+quark_string[i]+" flavour",samples,-0.5,hlimit,samples,-0.5,hlimit);
    for(int j=0;j<samples;j++){
      for(int k=0;k<samples;k++){	
	KS_d0bprob2[i]->SetBinContent(j+1,k+1,ptosigma(ks_d0bprob2[i][j][k]));
        KS_dEdxRatioPionOverKaonPri[i]->SetBinContent(j+1,k+1,ptosigma(ks_dEdxRatioPionOverKaonPri[i][j][k]));
	KS_dEdxRatioPionOverKaonSec[i]->SetBinContent(j+1,k+1,ptosigma(ks_dEdxRatioPionOverKaonSec[i][j][k]));
	KS_dEdxNKaonPri[i]->SetBinContent(j+1,k+1,ptosigma(ks_dEdxNKaonPri[i][j][k]));
	KS_dEdxNKaonSec[i]->SetBinContent(j+1,k+1,ptosigma(ks_dEdxNKaonSec[i][j][k]));
      }
    }
        
    for(int j=0;j<samples;j++){
      for(int k=0;k<samples;k++){
        Float_t d0b_b1=KS_d0bprob2[i]->GetBinContent(j+1,k+1);
        Float_t d0b_b2=KS_d0bprob2[i]->GetBinContent(k+1,j+1);
	Float_t d0b_bin=(d0b_b1+d0b_b2)/2;
	Float_t KoP_p_b1=KS_dEdxRatioPionOverKaonPri[i]->GetBinContent(j+1,k+1);
        Float_t KoP_p_b2=KS_dEdxRatioPionOverKaonPri[i]->GetBinContent(k+1,j+1);
        Float_t KoP_p_bin=(KoP_p_b1+KoP_p_b2)/2;
	Float_t KoP_s_b1=KS_dEdxRatioPionOverKaonSec[i]->GetBinContent(j+1,k+1);
	Float_t KoP_s_b2=KS_dEdxRatioPionOverKaonSec[i]->GetBinContent(k+1,j+1);
        Float_t KoP_s_bin=(KoP_s_b1+KoP_s_b2)/2;
	Float_t NK_p_b1=KS_dEdxNKaonPri[i]->GetBinContent(j+1,k+1);
	Float_t NK_p_b2=KS_dEdxNKaonPri[i]->GetBinContent(k+1,j+1);
        Float_t NK_p_bin=(NK_p_b1+NK_p_b2)/2;
	Float_t NK_s_b1=KS_dEdxNKaonSec[i]->GetBinContent(j+1,k+1);
        Float_t NK_s_b2=KS_dEdxNKaonSec[i]->GetBinContent(k+1,j+1);
        Float_t NK_s_bin=(NK_s_b1+NK_s_b2)/2;

	KS_d0bprob2[i]->SetBinContent(j+1,k+1,d0b_bin);
	KS_d0bprob2[i]->SetBinContent(k+1,j+1,d0b_bin);	
	KS_dEdxRatioPionOverKaonPri[i]->SetBinContent(j+1,k+1,KoP_p_bin);
        KS_dEdxRatioPionOverKaonPri[i]->SetBinContent(k+1,j+1,KoP_p_bin);
	KS_dEdxRatioPionOverKaonSec[i]->SetBinContent(j+1,k+1,KoP_s_bin);
	KS_dEdxRatioPionOverKaonSec[i]->SetBinContent(k+1,j+1,KoP_s_bin);
	KS_dEdxNKaonPri[i]->SetBinContent(j+1,k+1,NK_p_bin);
        KS_dEdxNKaonPri[i]->SetBinContent(k+1,j+1,NK_p_bin);
	KS_dEdxNKaonSec[i]->SetBinContent(j+1,k+1,NK_s_bin);
        KS_dEdxNKaonSec[i]->SetBinContent(k+1,j+1,NK_s_bin);
      }
    }
    
    KS_d0bprob2[i]->SetMaximum(5);
    KS_dEdxRatioPionOverKaonPri[i]->SetMaximum(5);
    KS_dEdxRatioPionOverKaonSec[i]->SetMaximum(5);
    KS_dEdxNKaonPri[i]->SetMaximum(5);
    KS_dEdxNKaonSec[i]->SetMaximum(5);
  }
  
  
  
  // PLOTS
  //SetQQbarStyle();
  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);
  //gStyle->SetOptTitle(0);
  //gStyle->SetTitleBorderSize(0);
  //gStyle->SetTitleStyle(0);
  //gStyle->SetTitleX(0.2);
  //gStyle->SetMarkerSize(1.5);
  

  //TCanvas *Kolmogorov_Test_ratios = new TCanvas("Kolmogorov_Test_ratios","Kolmogorov_Test_ratios",1200,600);
  //Kolmogorov_Test_ratios->Divide(3,1);
  TCanvas *Kolmogorov_Test_ratios = new TCanvas("Kolmogorov_Test_ratios","Kolmogorov_Test_ratios",900,900);
  Kolmogorov_Test_ratios->Divide(3,3);
  Kolmogorov_Test_ratios->cd(1);
  KS_d0bprob2[0]->SetTitle("KS Test (#sigma) for d0bprob2 (b-jets)");
  KS_d0bprob2[0]->GetYaxis()->SetTitle("Sample ID");
  KS_d0bprob2[0]->GetXaxis()->SetTitle("Sample ID");
  KS_d0bprob2[0]->Draw("COLZ");
  Kolmogorov_Test_ratios->cd(2);
  KS_d0bprob2[1]->SetTitle("KS Test (#sigma) for d0bprob2 (c-jets)");
  KS_d0bprob2[1]->GetYaxis()->SetTitle("Sample ID");
  KS_d0bprob2[1]->GetXaxis()->SetTitle("Sample ID");
  KS_d0bprob2[1]->Draw("COLZ");
  Kolmogorov_Test_ratios->cd(3);
  KS_d0bprob2[2]->SetTitle("KS Test (#sigma) for d0bprob2 (uds-jets)");
  KS_d0bprob2[2]->GetYaxis()->SetTitle("Sample ID");
  KS_d0bprob2[2]->GetXaxis()->SetTitle("Sample ID");
  KS_d0bprob2[2]->Draw("COLZ");
  
  //TCanvas *Kolmogorov_Test_dEdxRatioPionOverKaonPri = new TCanvas("Kolmogorov_Test_dEdxRatioPionOverKaonPri","Kolmogorov_Test_dEdxRatioPionOverKaonPri",1200,600);
  //Kolmogorov_Test_dEdxRatioPionOverKaonPri->Divide(3,1);
  //Kolmogorov_Test_dEdxRatioPionOverKaonPri->cd(1);
  Kolmogorov_Test_ratios->cd(4);
  KS_dEdxRatioPionOverKaonPri[0]->SetTitle("KS Test (#sigma) for dEdxRatioPionOverKaonPri (b-jets)");
  KS_dEdxRatioPionOverKaonPri[0]->GetYaxis()->SetTitle("Sample ID");
  KS_dEdxRatioPionOverKaonPri[0]->GetXaxis()->SetTitle("Sample ID");
  KS_dEdxRatioPionOverKaonPri[0]->Draw("COLZ");
  //Kolmogorov_Test_dEdxRatioPionOverKaonPri->cd(2);
  Kolmogorov_Test_ratios->cd(5);
  KS_dEdxRatioPionOverKaonPri[1]->SetTitle("KS Test (#sigma) for dEdxRatioPionOverKaonPri (c-jets)");
  KS_dEdxRatioPionOverKaonPri[1]->GetYaxis()->SetTitle("Sample ID");
  KS_dEdxRatioPionOverKaonPri[1]->GetXaxis()->SetTitle("Sample ID");
  KS_dEdxRatioPionOverKaonPri[1]->Draw("COLZ");
  //Kolmogorov_Test_dEdxRatioPionOverKaonPri->cd(3);
  Kolmogorov_Test_ratios->cd(6);
  KS_dEdxRatioPionOverKaonPri[2]->SetTitle("KS Test (#sigma) for dEdxRatioPionOverKaonPri (uds-jets)");
  KS_dEdxRatioPionOverKaonPri[2]->GetYaxis()->SetTitle("Sample ID");
  KS_dEdxRatioPionOverKaonPri[2]->GetXaxis()->SetTitle("Sample ID");
  KS_dEdxRatioPionOverKaonPri[2]->Draw("COLZ");

  //TCanvas *Kolmogorov_Test_dEdxRatioPionOverKaonSec = new TCanvas("Kolmogorov_Test_dEdxRatioPionOverKaonSec","Kolmogorov_Test_dEdxRatioPionOverKaonSec",1200,600);
  //Kolmogorov_Test_dEdxRatioPionOverKaonSec->Divide(3,1);
  //Kolmogorov_Test_dEdxRatioPionOverKaonSec->cd(1);
  Kolmogorov_Test_ratios->cd(7);
  KS_dEdxRatioPionOverKaonSec[0]->SetTitle("KS Test (#sigma) for dEdxRatioPionOverKaonSec (b-jets)");
  KS_dEdxRatioPionOverKaonSec[0]->GetYaxis()->SetTitle("Sample ID");
  KS_dEdxRatioPionOverKaonSec[0]->GetXaxis()->SetTitle("Sample ID");
  KS_dEdxRatioPionOverKaonSec[0]->Draw("COLZ");
  Kolmogorov_Test_ratios->cd(8);
  //Kolmogorov_Test_dEdxRatioPionOverKaonSec->cd(2);
  KS_dEdxRatioPionOverKaonSec[1]->SetTitle("KS Test (#sigma) for dEdxRatioPionOverKaonSec (c-jets)");
  KS_dEdxRatioPionOverKaonSec[1]->GetYaxis()->SetTitle("Sample ID");
  KS_dEdxRatioPionOverKaonSec[1]->GetXaxis()->SetTitle("Sample ID");
  KS_dEdxRatioPionOverKaonSec[1]->Draw("COLZ");
  //Kolmogorov_Test_dEdxRatioPionOverKaonSec->cd(3);
  Kolmogorov_Test_ratios->cd(9);
  KS_dEdxRatioPionOverKaonSec[2]->SetTitle("KS Test (#sigma) for dEdxRatioPionOverKaonSec (uds-jets)");
  KS_dEdxRatioPionOverKaonSec[2]->GetYaxis()->SetTitle("Sample ID");
  KS_dEdxRatioPionOverKaonSec[2]->GetXaxis()->SetTitle("Sample ID");
  KS_dEdxRatioPionOverKaonSec[2]->Draw("COLZ");
  
  
  TCanvas *Kolmogorov_Test_nkaon = new TCanvas("Kolmogorov_Test_nkaon","Kolmogorov_Test_nkaon",900,900);
  Kolmogorov_Test_nkaon->Divide(3,3);
  Kolmogorov_Test_nkaon->cd(1);
  KS_d0bprob2[0]->SetTitle("KS Test (#sigma) for d0bprob2 (b-jets)");
  KS_d0bprob2[0]->GetYaxis()->SetTitle("Sample ID");
  KS_d0bprob2[0]->GetXaxis()->SetTitle("Sample ID");
  KS_d0bprob2[0]->Draw("COLZ");
  Kolmogorov_Test_nkaon->cd(2);
  KS_d0bprob2[1]->SetTitle("KS Test (#sigma) for d0bprob2 (c-jets)");
  KS_d0bprob2[1]->GetYaxis()->SetTitle("Sample ID");
  KS_d0bprob2[1]->GetXaxis()->SetTitle("Sample ID");
  KS_d0bprob2[1]->Draw("COLZ");
  Kolmogorov_Test_nkaon->cd(3);
  KS_d0bprob2[2]->SetTitle("KS Test (#sigma) for d0bprob2 (uds-jets)");
  KS_d0bprob2[2]->GetYaxis()->SetTitle("Sample ID");
  KS_d0bprob2[2]->GetXaxis()->SetTitle("Sample ID");
  KS_d0bprob2[2]->Draw("COLZ");
  Kolmogorov_Test_nkaon->cd(4);
  KS_dEdxNKaonPri[0]->SetTitle("KS Test (#sigma) for dEdxNKaonPri (b-jets)");
  KS_dEdxNKaonPri[0]->GetYaxis()->SetTitle("Sample ID");
  KS_dEdxNKaonPri[0]->GetXaxis()->SetTitle("Sample ID");
  KS_dEdxNKaonPri[0]->Draw("COLZ");
  Kolmogorov_Test_nkaon->cd(5);
  KS_dEdxNKaonPri[1]->SetTitle("KS Test (#sigma) for dEdxNKaonPri (c-jets)");
  KS_dEdxNKaonPri[1]->GetYaxis()->SetTitle("Sample ID");
  KS_dEdxNKaonPri[1]->GetXaxis()->SetTitle("Sample ID");
  KS_dEdxNKaonPri[1]->Draw("COLZ");
  Kolmogorov_Test_nkaon->cd(6);
  KS_dEdxNKaonPri[2]->SetTitle("KS Test (#sigma) for dEdxNKaonPri (uds-jets)");
  KS_dEdxNKaonPri[2]->GetYaxis()->SetTitle("Sample ID");
  KS_dEdxNKaonPri[2]->GetXaxis()->SetTitle("Sample ID");
  KS_dEdxNKaonPri[2]->Draw("COLZ");
  Kolmogorov_Test_nkaon->cd(7);
  KS_dEdxNKaonSec[0]->SetTitle("KS Test (#sigma) for dEdxNKaonSec (b-jets)");
  KS_dEdxNKaonSec[0]->GetYaxis()->SetTitle("Sample ID");
  KS_dEdxNKaonSec[0]->GetXaxis()->SetTitle("Sample ID");
  KS_dEdxNKaonSec[0]->Draw("COLZ");
  Kolmogorov_Test_nkaon->cd(8);
  KS_dEdxNKaonSec[1]->SetTitle("KS Test (#sigma) for dEdxNKaonSec (c-jets)");
  KS_dEdxNKaonSec[1]->GetYaxis()->SetTitle("Sample ID");
  KS_dEdxNKaonSec[1]->GetXaxis()->SetTitle("Sample ID");
  KS_dEdxNKaonSec[1]->Draw("COLZ");
  Kolmogorov_Test_nkaon->cd(9);
  KS_dEdxNKaonSec[2]->SetTitle("KS Test (#sigma) for dEdxNKaonSec (uds-jets)");
  KS_dEdxNKaonSec[2]->GetYaxis()->SetTitle("Sample ID");
  KS_dEdxNKaonSec[2]->GetXaxis()->SetTitle("Sample ID");
  KS_dEdxNKaonSec[2]->Draw("COLZ");
}
