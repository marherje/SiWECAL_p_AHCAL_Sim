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
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystemFile.h"
#include "style/Style.C"
#include "style/Labels.C"

void lcfiplusvariables(TString energy, TString pol, TString cat){
  
  // Data:
  std::cout<<"Init Reading"<<std::endl;
  std::cout<<energy<<endl;
  
  TString lumi;
  if(energy=="250"){
    if(pol=="eL_pR")lumi="20";
    if(pol=="eR_pL")lumi="36";
  }
  if(energy=="500"){
    if(pol=="eL_pR")lumi="1430";
    if(pol=="eR_pL")lumi="1060";
    //    if(pol=="eL_pR")lumi="47";
    //if(pol=="eR_pL")lumi="47";
  }
  
  //folder: \"MakeNTuples500/\"
  TString foldername="../../../FlavourTagExtension/ReTraining500NewSamples/MakeNTuples500_beta/";
  //TString foldername="../../data/MakeNTuples"+energy+"/";
  TString quark_string[3]={"b_quark","c_quark","light_quark"};
  TString filename[3];
  TFile *f[3];
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
  TBranch *trk1d0sig;
  TBranch *trk2d0sig;
  TBranch *trk1z0sig;
  TBranch *trk2z0sig;
  TBranch *trk1pt;
  TBranch *trk2pt;
  TBranch *jprobr;
  TBranch *jprobr5sigma;
  TBranch *jprobz;
  TBranch *jprobz5sigma;
  TBranch *d0bprob2;
  TBranch *d0cprob2;
  TBranch *d0qprob2;
  TBranch *z0bprob2;
  TBranch *z0cprob2;
  TBranch *z0qprob2;
  TBranch *nmuon;
  TBranch *nelectron;
  TBranch *trkmass;
  TBranch *dEdxRatioPionOverKaonPri;
  TBranch *dEdxRatioPionOverKaonSec;
  //  if(cat == "B" || cat == "C" || cat == "D"){
    TBranch *one_vtxprob; //BCD
    TBranch *vtxlen1; //BCD
    TBranch *vtxsig1; //BCD
    TBranch *vtxdirang1; //BCD
    TBranch *vtxmult1; //BCD
    TBranch *vtxmom1; //BCD
    TBranch *vtxmass1; //BCD
    TBranch *vtxmass; //BCD
    TBranch *vtxmasspc; //BCD
    TBranch *vtxprob; //BCD
    //}
    //if(cat=="D"){
    TBranch *vtxlen2; //D
    TBranch *vtxlen12; //D
    TBranch *vtxsig2; //D
    TBranch *vtxsig12; //D
    TBranch *vtxdirang2; //D
    TBranch *vtxmult2; //D
    TBranch *vtxmult; //D
    TBranch *vtxmom2; //D
    TBranch *vtxmass2; //D
    //}
  // Addresses
  Float_t nvtx_a;
  Float_t nvtxall_a; //control
  Float_t trk1d0sig_a;
  Float_t trk2d0sig_a;
  Float_t trk1z0sig_a;
  Float_t trk2z0sig_a;
  Float_t trk1pt_a;
  Float_t trk2pt_a;
  Float_t jprobr_a;
  Float_t jprobr5sigma_a;
  Float_t jprobz_a;
  Float_t jprobz5sigma_a;
  Float_t d0bprob2_a;
  Float_t d0cprob2_a;
  Float_t d0qprob2_a;
  Float_t z0bprob2_a;
  Float_t z0cprob2_a;
  Float_t z0qprob2_a;
  Float_t nmuon_a;
  Float_t nelectron_a;
  Float_t trkmass_a;
  Float_t dEdxRatioPionOverKaonPri_a;
  Float_t dEdxRatioPionOverKaonSec_a;
  //if(cat == "B" || cat == "C" || cat == "D"){
    Float_t one_vtxprob_a; //BCD
    Float_t vtxlen1_a; //BCD
    Float_t vtxsig1_a; //BCD
    Float_t vtxdirang1_a; //BCD
    Float_t vtxmult1_a; //BCD
    Float_t vtxmom1_a; //BCD
    Float_t vtxmass1_a; //BCD
    Float_t vtxmass_a; //BCD
    Float_t vtxmasspc_a; //BCD
    Float_t vtxprob_a; //BCD
    //}
    //if(cat=="D"){
    Float_t vtxlen2_a; //D
    Float_t vtxlen12_a; //D
    Float_t vtxsig2_a; //D
    Float_t vtxsig12_a; //D
    Float_t vtxdirang2_a; //D
    Float_t vtxmult2_a; //D
    Float_t vtxmult_a; //D
    Float_t vtxmom2_a; //D
    Float_t vtxmass2_a; //D
    //}

  // To store and plot:
  TH1F *trk1d0sig_histo[3];
  TCanvas *trk1d0sig_canvas = new TCanvas("trk1d0sig canvas","trk1d0sig canvas",1000,1000);
  TH1F *trk2d0sig_histo[3];
  TCanvas *trk2d0sig_canvas = new TCanvas("trk2d0sig canvas","trk2d0sig canvas",1000,1000);
  TH1F *trk1z0sig_histo[3];
  TCanvas *trk1z0sig_canvas = new TCanvas("trk1z0sig canvas","trk1z0sig canvas",1000,1000);
  TH1F *trk2z0sig_histo[3];
  TCanvas *trk2z0sig_canvas = new TCanvas("trk2z0sig canvas","trk2z0sig canvas",1000,1000);
  TH1F *trk1pt_histo[3];
  TCanvas *trk1pt_canvas = new TCanvas("trk1pt canvas","trk1pt canvas",1000,1000);
  TH1F *trk2pt_histo[3];
  TCanvas *trk2pt_canvas = new TCanvas("trk2pt canvas","trk2pt canvas",1000,1000);
  TH1F *jprobr_histo[3];
  TCanvas *jprobr_canvas = new TCanvas("jprobr canvas","jprobr canvas",1000,1000);
  TH1F *jprobr5sigma_histo[3];
  TCanvas *jprobr5sigma_canvas = new TCanvas("jprobr5sigma canvas","jprobr5sigma canvas",1000,1000);
  TH1F *jprobz_histo[3];
  TCanvas *jprobz_canvas = new TCanvas("jprobz canvas","jprobz canvas",1000,1000);
  TH1F *jprobz5sigma_histo[3];
  TCanvas *jprobz5sigma_canvas = new TCanvas("jprobz5sigma canvas","jprobz5sigma canvas",1000,1000);
  TH1F *d0bprob2_histo[3];
  TCanvas *d0bprob2_canvas = new TCanvas("d0bprob2 canvas","d0bprob2 canvas",1000,1000);
  TH1F *d0cprob2_histo[3];
  TCanvas *d0cprob2_canvas = new TCanvas("d0cprob2 canvas","d0cprob2 canvas",1000,1000);
  TH1F *d0qprob2_histo[3];
  TCanvas *d0qprob2_canvas = new TCanvas("d0qprob2 canvas","d0qprob2 canvas",1000,1000);
  TH1F *z0bprob2_histo[3];
  TCanvas *z0bprob2_canvas = new TCanvas("z0bprob2 canvas","z0bprob2 canvas",1000,1000);
  TH1F *z0cprob2_histo[3];
  TCanvas *z0cprob2_canvas = new TCanvas("z0cprob2 canvas","z0cprob2 canvas",1000,1000);
  TH1F *z0qprob2_histo[3];
  TCanvas *z0qprob2_canvas = new TCanvas("z0qprob2 canvas","z0qprob2 canvas",1000,1000);
  TH1F *nmuon_histo[3];
  TCanvas *nmuon_canvas = new TCanvas("nmuon canvas","nmuon canvas",1000,1000);
  TH1F *nelectron_histo[3];
  TCanvas *nelectron_canvas = new TCanvas("nelectron canvas","nelectron canvas",1000,1000);
  TH1F *trkmass_histo[3];
  TCanvas *trkmass_canvas = new TCanvas("trkmass canvas","trkmass canvas",1000,1000);
  TH1F *dEdxRatioPionOverKaonPri_histo[3];
  TCanvas *dEdxRatioPionOverKaonPri_canvas = new TCanvas("dEdxRatioPionOverKaonPri canvas","dEdxRatioPionOverKaonPri canvas",1000,1000);
  TH1F *dEdxRatioPionOverKaonSec_histo[3];
  TCanvas *dEdxRatioPionOverKaonSec_canvas = new TCanvas("dEdxRatioPionOverKaonSec canvas","dEdxRatioPionOverKaonSec canvas",1000,1000);
  
  //BCD
  //  if(cat == "B" || cat == "C" || cat == "D"){
    TH1F *one_vtxprob_histo[3];
    TCanvas *one_vtxprob_canvas = new TCanvas("one_vtxprob canvas","one_vtxprob canvas",1000,1000);
    TH1F *vtxlen1_histo[3];
    TCanvas *vtxlen1_canvas = new TCanvas("vtxlen1 canvas","vtxlen1 canvas",1000,1000);
    TH1F *vtxsig1_histo[3];
    TCanvas *vtxsig1_canvas = new TCanvas("vtxsig1 canvas","vtxsig1 canvas",1000,1000);
    TH1F *vtxdirang1_histo[3];
    TCanvas *vtxdirang1_canvas = new TCanvas("vtxdirang1 canvas","vtxdirang1 canvas",1000,1000);
    TH1F *vtxmult1_histo[3];
    TCanvas *vtxmult1_canvas = new TCanvas("vtxmult1 canvas","vtxmult1 canvas",1000,1000);
    TH1F *vtxmom1_histo[3];
    TCanvas *vtxmom1_canvas = new TCanvas("vtxmom1 canvas","vtxmom1 canvas",1000,1000);
    TH1F *vtxmass1_histo[3];
    TCanvas *vtxmass1_canvas = new TCanvas("vtxmass1 canvas","vtxmass1 canvas",1000,1000);
    TH1F *vtxmass_histo[3];
    TCanvas *vtxmass_canvas = new TCanvas("vtxmass canvas","vtxmass canvas",1000,1000);
    TH1F *vtxmasspc_histo[3];
    TCanvas *vtxmasspc_canvas = new TCanvas("vtxmasspc canvas","vtxmasspc canvas",1000,1000);
    TH1F *vtxprob_histo[3];
    TCanvas *vtxprob_canvas = new TCanvas("vtxprob canvas","vtxprob canvas",1000,1000);
    //}
  //D
    //  if(cat=="D"){
    TH1F *vtxlen2_histo[3];
    TCanvas *vtxlen2_canvas = new TCanvas("vtxlen2 canvas","vtxlen2 canvas",1000,1000);
    TH1F *vtxlen12_histo[3];
    TCanvas *vtxlen12_canvas = new TCanvas("vtxlen12 canvas","vtxlen12 canvas",1000,1000);
    TH1F *vtxsig2_histo[3];
    TCanvas *vtxsig2_canvas = new TCanvas("vtxsig2 canvas","vtxsig2 canvas",1000,1000);
    TH1F *vtxsig12_histo[3];
    TCanvas *vtxsig12_canvas = new TCanvas("vtxsig12 canvas","vtxsig12 canvas",1000,1000);
    TH1F *vtxdirang2_histo[3];
    TCanvas *vtxdirang2_canvas = new TCanvas("vtxdirang2 canvas","vtxdirang2 canvas",1000,1000);
    TH1F *vtxmult2_histo[3];
    TCanvas *vtxmult2_canvas = new TCanvas("vtxmult2 canvas","vtxmult2 canvas",1000,1000);
    TH1F *vtxmult_histo[3];
    TCanvas *vtxmult_canvas = new TCanvas("vtxmult canvas","vtxmult canvas",1000,1000);
    TH1F *vtxmom2_histo[3];
    TCanvas *vtxmom2_canvas = new TCanvas("vtxmom2 canvas","vtxmom2 canvas",1000,1000);
    TH1F *vtxmass2_histo[3];
    TCanvas *vtxmass2_canvas = new TCanvas("vtxmass2 canvas","vtxmass2 canvas",1000,1000);
    //}

  for(int i=0;i<3;i++){
    filename[i]=foldername+pol+"_"+quark_string[i]+"_merged.root";
    //eL_pR_b_quark_merged.root
    f[i]=new TFile(filename[i]);
    std::cout<<filename[i]<<std::endl;

    // Load tree
    tree = f[i]->Get<TTree>("ntp");
    // Load branches
    nvtx = tree->GetBranch("nvtx");
    nvtxall = tree->GetBranch("nvtxall");
    
    trk1d0sig = tree->GetBranch("trk1d0sig");
    trk2d0sig = tree->GetBranch("trk2d0sig");
    trk1z0sig = tree->GetBranch("trk1z0sig");
    trk2z0sig = tree->GetBranch("trk2z0sig");
    trk1pt = tree->GetBranch("trk1pt");
    trk2pt = tree->GetBranch("trk2pt");
    jprobr = tree->GetBranch("jprobr");
    jprobr5sigma = tree->GetBranch("jprobr5sigma");
    jprobz = tree->GetBranch("jprobz");
    jprobz5sigma = tree->GetBranch("jprobz5sigma");
    d0bprob2 = tree->GetBranch("d0bprob2");
    d0cprob2 = tree->GetBranch("d0cprob2");
    d0qprob2 = tree->GetBranch("d0qprob2");
    z0bprob2 = tree->GetBranch("z0bprob2");
    z0cprob2 = tree->GetBranch("z0cprob2");
    z0qprob2 = tree->GetBranch("z0qprob2");
    nmuon = tree->GetBranch("nmuon");
    nelectron = tree->GetBranch("nelectron");
    trkmass = tree->GetBranch("trkmass");
    dEdxRatioPionOverKaonPri = tree->GetBranch("dEdxRatioPionOverKaonPri");
    dEdxRatioPionOverKaonSec = tree->GetBranch("dEdxRatioPionOverKaonSec");

    //BCD
    if(cat == "B" || cat == "C" || cat == "D"){
      one_vtxprob = tree->GetBranch("1vtxprob");
      vtxlen1 = tree->GetBranch("vtxlen1");
      vtxsig1 = tree->GetBranch("vtxsig1");
      vtxdirang1 = tree->GetBranch("vtxdirang1");
      vtxmult1 = tree->GetBranch("vtxmult1");
      vtxmom1 = tree->GetBranch("vtxmom1");
      vtxmass1 = tree->GetBranch("vtxmass1");
      vtxmass = tree->GetBranch("vtxmass");
      vtxmasspc = tree->GetBranch("vtxmasspc");
      vtxprob = tree->GetBranch("vtxprob");
    }
    //D
    if(cat=="D"){
      vtxlen2 = tree->GetBranch("vtxlen2");
      vtxlen12 = tree->GetBranch("vtxlen12");
      vtxsig2 = tree->GetBranch("vtxsig2");
      vtxsig12 = tree->GetBranch("vtxsig12");
      vtxdirang2 = tree->GetBranch("vtxdirang2");
      vtxmult2 = tree->GetBranch("vtxmult2");
      vtxmult = tree->GetBranch("vtxmult");
      vtxmom2 = tree->GetBranch("vtxmom2");
      vtxmass2 = tree->GetBranch("vtxmass2");
    }
    // Set Addresses
    nvtx->SetAddress(&nvtx_a);
    nvtxall->SetAddress(&nvtxall_a);
    
    trk1d0sig->SetAddress(&trk1d0sig_a);
    trk2d0sig->SetAddress(&trk2d0sig_a);
    trk1z0sig->SetAddress(&trk1z0sig_a);
    trk2z0sig->SetAddress(&trk2z0sig_a);
    trk1pt->SetAddress(&trk1pt_a);
    trk2pt->SetAddress(&trk2pt_a);
    jprobr->SetAddress(&jprobr_a);
    jprobr5sigma->SetAddress(&jprobr5sigma_a);
    jprobz->SetAddress(&jprobz_a);
    jprobz5sigma->SetAddress(&jprobz5sigma_a);
    d0bprob2->SetAddress(&d0bprob2_a);
    d0cprob2->SetAddress(&d0cprob2_a);
    d0qprob2->SetAddress(&d0qprob2_a);
    z0bprob2->SetAddress(&z0bprob2_a);
    z0cprob2->SetAddress(&z0cprob2_a);
    z0qprob2->SetAddress(&z0qprob2_a);
    nmuon->SetAddress(&nmuon_a);
    nelectron->SetAddress(&nelectron_a);
    trkmass->SetAddress(&trkmass_a);
    dEdxRatioPionOverKaonPri->SetAddress(&dEdxRatioPionOverKaonPri_a);
    dEdxRatioPionOverKaonSec->SetAddress(&dEdxRatioPionOverKaonSec_a);

    //BCD
    if(cat!="A"){
      one_vtxprob->SetAddress(&one_vtxprob_a);
      vtxlen1->SetAddress(&vtxlen1_a);
      vtxsig1->SetAddress(&vtxsig1_a);
      vtxdirang1->SetAddress(&vtxdirang1_a);
      vtxmult1->SetAddress(&vtxmult1_a);
      vtxmom1->SetAddress(&vtxmom1_a);
      vtxmass1->SetAddress(&vtxmass1_a);
      vtxmass->SetAddress(&vtxmass_a);
      vtxmasspc->SetAddress(&vtxmasspc_a);
      vtxprob->SetAddress(&vtxprob_a);
    }
    //D
    if(cat=="D"){
      vtxlen2->SetAddress(&vtxlen2_a);
      vtxlen12->SetAddress(&vtxlen12_a);
      vtxsig2->SetAddress(&vtxsig2_a);
      vtxsig12->SetAddress(&vtxsig12_a);
      vtxdirang2->SetAddress(&vtxdirang2_a);
      vtxmult2->SetAddress(&vtxmult2_a);
      vtxmult->SetAddress(&vtxmult_a);
      vtxmom2->SetAddress(&vtxmom2_a);
      vtxmass2->SetAddress(&vtxmass2_a);
    }
    int totalentries;
    totalentries = d0bprob2->GetEntries();    

    Float_t xmin_trk1d0sig = tree->GetMinimum("trk1d0sig");
    Float_t xmax_trk1d0sig = tree->GetMaximum("trk1d0sig");
    Float_t xmin_trk2d0sig = tree->GetMinimum("trk2d0sig");
    Float_t xmax_trk2d0sig = tree->GetMaximum("trk2d0sig");
    Float_t xmin_trk1z0sig = tree->GetMinimum("trk1z0sig");
    Float_t xmax_trk1z0sig = tree->GetMaximum("trk1z0sig");
    Float_t xmin_trk2z0sig = tree->GetMinimum("trk2z0sig");
    Float_t xmax_trk2z0sig = tree->GetMaximum("trk2z0sig");
    Float_t xmin_trk1pt = tree->GetMinimum("trk1pt");
    Float_t xmax_trk1pt = tree->GetMaximum("trk1pt");
    Float_t xmin_trk2pt = tree->GetMinimum("trk2pt");
    Float_t xmax_trk2pt = tree->GetMaximum("trk2pt");
    Float_t xmin_jprobr = tree->GetMinimum("jprobr");
    Float_t xmax_jprobr = tree->GetMaximum("jprobr");
    Float_t xmin_jprobr5sigma = tree->GetMinimum("jprobr5sigma");
    Float_t xmax_jprobr5sigma = tree->GetMaximum("jprobr5sigma");
    Float_t xmin_jprobz = tree->GetMinimum("jprobz");
    Float_t xmax_jprobz = tree->GetMaximum("jprobz");
    Float_t xmin_jprobz5sigma = tree->GetMinimum("jprobz5sigma");
    Float_t xmax_jprobz5sigma = tree->GetMaximum("jprobz5sigma");

    Float_t xmin_d0bprob2 = tree->GetMinimum("d0bprob2");
    Float_t xmax_d0bprob2 = tree->GetMaximum("d0bprob2");
    Float_t xmin_d0cprob2 = tree->GetMinimum("d0cprob2");
    Float_t xmax_d0cprob2 = tree->GetMaximum("d0cprob2");
    Float_t xmin_d0qprob2 = tree->GetMinimum("d0qprob2");
    Float_t xmax_d0qprob2 = tree->GetMaximum("d0qprob2");
    Float_t xmin_z0bprob2 = tree->GetMinimum("z0bprob2");
    Float_t xmax_z0bprob2 = tree->GetMaximum("z0bprob2");
    Float_t xmin_z0cprob2 = tree->GetMinimum("z0cprob2");
    Float_t xmax_z0cprob2 = tree->GetMaximum("z0cprob2");
    Float_t xmin_z0qprob2 = tree->GetMinimum("z0qprob2");
    Float_t xmax_z0qprob2 = tree->GetMaximum("z0qprob2");
    Float_t xmin_nmuon = tree->GetMinimum("nmuon");
    Float_t xmax_nmuon = tree->GetMaximum("nmuon");
    Float_t xmin_nelectron = tree->GetMinimum("nelectron");
    Float_t xmax_nelectron = tree->GetMaximum("nelectron");
    Float_t xmin_trkmass = tree->GetMinimum("trkmass");
    Float_t xmax_trkmass = tree->GetMaximum("trkmass");
    Float_t xmin_dEdxRatioPionOverKaonPri= tree->GetMinimum("dEdxRatioPionOverKaonPri");
    Float_t xmax_dEdxRatioPionOverKaonPri= tree->GetMaximum("dEdxRatioPionOverKaonPri");
    Float_t xmin_dEdxRatioPionOverKaonSec= tree->GetMinimum("dEdxRatioPionOverKaonSec");
    Float_t xmax_dEdxRatioPionOverKaonSec= tree->GetMaximum("dEdxRatioPionOverKaonSec");

    //BCD
    if(cat!="A"){
      Float_t xmin_one_vtxprob = tree->GetMinimum("one_vtxprob");
      Float_t xmax_one_vtxprob = tree->GetMaximum("one_vtxprob");
      Float_t xmin_vtxlen1 = tree->GetMinimum("vtxlen1");
      Float_t xmax_vtxlen1 = tree->GetMaximum("vtxlen1");
      Float_t xmin_vtxsig1 = tree->GetMinimum("vtxsig1");
      Float_t xmax_vtxsig1 = tree->GetMaximum("vtxsig1");
      Float_t xmin_vtxdirang1 = tree->GetMinimum("vtxdirang1");
      Float_t xmax_vtxdirang1 = tree->GetMaximum("vtxdirang1");
      Float_t xmin_vtxmult1 = tree->GetMinimum("vtxmult1");
      Float_t xmax_vtxmult1 = tree->GetMaximum("vtxmult1");
      Float_t xmin_vtxmom1 = tree->GetMinimum("vtxmom1");
      Float_t xmax_vtxmom1 = tree->GetMaximum("vtxmom1");
      Float_t xmin_vtxmass1 = tree->GetMinimum("vtxmass1");
      Float_t xmax_vtxmass1 = tree->GetMaximum("vtxmass1");
      Float_t xmin_vtxmass = tree->GetMinimum("vtxmass");
      Float_t xmax_vtxmass = tree->GetMaximum("vtxmass");
      Float_t xmin_vtxmasspc = tree->GetMinimum("vtxmasspc");
      Float_t xmax_vtxmasspc = tree->GetMaximum("vtxmasspc");
      Float_t xmin_vtxprob = tree->GetMinimum("vtxprob");
      Float_t xmax_vtxprob = tree->GetMaximum("vtxprob");
    }
    //D
    if(cat=="D"){
      Float_t xmin_vtxlen2 = tree->GetMinimum("vtxlen2");
      Float_t xmax_vtxlen2 = tree->GetMaximum("vtxlen2");
      Float_t xmin_vtxlen12 = tree->GetMinimum("vtxlen12");
      Float_t xmax_vtxlen12 = tree->GetMaximum("vtxlen12");
      Float_t xmin_vtxsig2 = tree->GetMinimum("vtxsig2");
      Float_t xmax_vtxsig2 = tree->GetMaximum("vtxsig2");
      Float_t xmin_vtxsig12 = tree->GetMinimum("vtxsig12");
      Float_t xmax_vtxsig12 = tree->GetMaximum("vtxsig12");
      Float_t xmin_vtxdirang2 = tree->GetMinimum("vtxdirang2");
      Float_t xmax_vtxdirang2 = tree->GetMaximum("vtxdirang2");
      Float_t xmin_vtxmult2 = tree->GetMinimum("vtxmult2");
      Float_t xmax_vtxmult2 = tree->GetMaximum("vtxmult2");
      Float_t xmin_vtxmult = tree->GetMinimum("vtxmult");
      Float_t xmax_vtxmult = tree->GetMaximum("vtxmult");
      Float_t xmin_vtxmom2 = tree->GetMinimum("vtxmom2");
      Float_t xmax_vtxmom2 = tree->GetMaximum("vtxmom2");
      Float_t xmin_vtxmass2 = tree->GetMinimum("vtxmass2");
      Float_t xmax_vtxmass2 = tree->GetMinimum("vtxmass2");    
    }
      std::cout<<"Entries in category "<<cat<<": "<<totalentries<<std::endl;
    //std::cout<<"Min: "<<xmin_d0bprob2<<std::endl;
    //std::cout<<"Max: "<<xmax_d0bprob2<<std::endl;
    
    trk1d0sig_histo[i] = new TH1F("trk1d0sig "+quark_string[i],"trk1d0sig "+quark_string[i],100,0,5000);
    trk2d0sig_histo[i] = new TH1F("trk2d0sig "+quark_string[i],"trk2d0sig "+quark_string[i],100,0,5000);
    trk1z0sig_histo[i] = new TH1F("trk1z0sig "+quark_string[i],"trk1z0sig "+quark_string[i],100,-50,50);
    trk2z0sig_histo[i] = new TH1F("trk2z0sig "+quark_string[i],"trk2z0sig "+quark_string[i],100,-500,500);
    trk1pt_histo[i] = new TH1F("trk1pt "+quark_string[i],"trk1pt "+quark_string[i],100,0,100);
    trk2pt_histo[i] = new TH1F("trk2pt "+quark_string[i],"trk2pt "+quark_string[i],100,0,100);
    jprobr_histo[i] = new TH1F("jprobr "+quark_string[i],"jprobr "+quark_string[i],100,0,1);
    jprobz_histo[i] = new TH1F("jprobz "+quark_string[i],"jprobz "+quark_string[i],100,0,1);
    jprobr5sigma_histo[i] = new TH1F("jprobr5sigma "+quark_string[i],"jprobr5sigma "+quark_string[i],100,0,1);
    jprobz5sigma_histo[i] = new TH1F("jprobz5sigma "+quark_string[i],"jprobz5sigma "+quark_string[i],100,0,1);
    d0bprob2_histo[i] = new TH1F("d0bprob2 "+quark_string[i],"d0bprob2 "+quark_string[i],100,0,1);
    d0cprob2_histo[i] = new TH1F("d0cprob2 "+quark_string[i],"d0cprob2 "+quark_string[i],100,0,1);
    d0qprob2_histo[i] = new TH1F("d0qprob2 "+quark_string[i],"d0qprob2 "+quark_string[i],100,0,1);
    z0bprob2_histo[i] = new TH1F("z0bprob2 "+quark_string[i],"z0bprob2 "+quark_string[i],100,0,1);
    z0cprob2_histo[i] = new TH1F("z0cprob2 "+quark_string[i],"z0cprob2 "+quark_string[i],100,0,1);
    z0qprob2_histo[i] = new TH1F("z0qprob2 "+quark_string[i],"z0qprob2 "+quark_string[i],100,0,1);
    nmuon_histo[i] = new TH1F("nmuon "+quark_string[i],"nmuon "+quark_string[i],100,0,3);
    nelectron_histo[i] = new TH1F("nelectron "+quark_string[i],"nelectron "+quark_string[i],100,0,4);
    trkmass_histo[i] = new TH1F("trkmass "+quark_string[i],"trkmass "+quark_string[i],100,0,20);
    dEdxRatioPionOverKaonPri_histo[i] = new TH1F("dEdxRatioPionOverKaonPri "+quark_string[i],"dEdxRatioPionOverKaonPri "+quark_string[i],100,-3,10);
    dEdxRatioPionOverKaonSec_histo[i] = new TH1F("dEdxRatioPionOverKaonSec "+quark_string[i],"dEdxRatioPionOverKaonSec "+quark_string[i],100,-3,10);
    //BCD
    if(cat!="A"){
      one_vtxprob_histo[i] = new TH1F("one_vtxprob "+quark_string[i],"one_vtxprob "+quark_string[i],100,0,0.00001);
      vtxlen1_histo[i] = new TH1F("vtxlen1 "+quark_string[i],"vtxlen1 "+quark_string[i],100,0,200);
      vtxsig1_histo[i] = new TH1F("vtxsig1 "+quark_string[i],"vtxsig1 "+quark_string[i],100,0,5000);
      vtxdirang1_histo[i] = new TH1F("vtxdirang1 "+quark_string[i],"vtxdirang1 "+quark_string[i],100,0,3.15);
      vtxmult1_histo[i] = new TH1F("vtxmult1 "+quark_string[i],"vtxmult1 "+quark_string[i],100,0,16);
      vtxmom1_histo[i] = new TH1F("vtxmom1 "+quark_string[i],"vtxmom1 "+quark_string[i],100,0,200);
      vtxmass1_histo[i] = new TH1F("vtxmass1 "+quark_string[i],"vtxmass1 "+quark_string[i],100,0,20);
      vtxmass_histo[i] = new TH1F("vtxmass "+quark_string[i],"vtxmass "+quark_string[i],100,0,20);
      vtxmasspc_histo[i] = new TH1F("vtxmasspc "+quark_string[i],"vtxmasspc "+quark_string[i],100,0,20);
      vtxprob_histo[i] = new TH1F("vtxprob "+quark_string[i],"vtxprob "+quark_string[i],100,0,1);
    }
    //D
    if(cat=="D"){
      vtxlen2_histo[i] = new TH1F("vtxlen2 "+quark_string[i],"vtxlen2 "+quark_string[i],100,0,200);
      vtxlen12_histo[i] = new TH1F("vtxlen12 "+quark_string[i],"vtxlen12 "+quark_string[i],100,0,35);
      vtxsig2_histo[i] = new TH1F("vtxsig2 "+quark_string[i],"vtxsig2 "+quark_string[i],100,0,8000);
      vtxsig12_histo[i] = new TH1F("vtxsig12 "+quark_string[i],"vtxsig12 "+quark_string[i],100,0,4000);
      vtxdirang2_histo[i] = new TH1F("vtxdirang2 "+quark_string[i],"vtxdirang2 "+quark_string[i],100,0,3.15);
      vtxmult2_histo[i] = new TH1F("vtxmult2 "+quark_string[i],"vtxmult2 "+quark_string[i],100,0,12);
      vtxmult_histo[i] = new TH1F("vtxmult "+quark_string[i],"vtxmult "+quark_string[i],100,0,22);
      vtxmom2_histo[i] = new TH1F("vtxmom2 "+quark_string[i],"vtxmom2 "+quark_string[i],100,0,200);
      vtxmass2_histo[i] = new TH1F("vtxmass2 "+quark_string[i],"vtxmass2 "+quark_string[i],100,0,10);
    }

    for(int ientry=0;ientry<totalentries;ientry++){
      nvtx->GetEntry(ientry);
      nvtxall->GetEntry(ientry);
      //if(nvtx_a==0)cout<<"vtx: "<<nvtx_a<<"Total vtx: "<<nvtxall_a<<endl;
      // Variables
      trk1d0sig->GetEntry(ientry);
      trk2d0sig->GetEntry(ientry);
      trk1z0sig->GetEntry(ientry);
      trk2z0sig->GetEntry(ientry);
      trk1pt->GetEntry(ientry);
      trk2pt->GetEntry(ientry);
      jprobr->GetEntry(ientry);
      jprobr5sigma->GetEntry(ientry);
      jprobz->GetEntry(ientry);
      jprobz5sigma->GetEntry(ientry);
      d0bprob2->GetEntry(ientry);
      d0cprob2->GetEntry(ientry);
      d0qprob2->GetEntry(ientry);
      z0bprob2->GetEntry(ientry);
      z0cprob2->GetEntry(ientry);
      z0qprob2->GetEntry(ientry);
      nmuon->GetEntry(ientry);
      nelectron->GetEntry(ientry);
      trkmass->GetEntry(ientry);
      dEdxRatioPionOverKaonPri->GetEntry(ientry);
      dEdxRatioPionOverKaonSec->GetEntry(ientry);
      //BCD
      if(cat == "B" || cat == "C" || cat == "D"){
	one_vtxprob->GetEntry(ientry);
	vtxlen1->GetEntry(ientry);
	vtxsig1->GetEntry(ientry);
	vtxdirang1->GetEntry(ientry);
	vtxmult1->GetEntry(ientry);
	vtxmom1->GetEntry(ientry);
	vtxmass1->GetEntry(ientry);
	vtxmass->GetEntry(ientry);
	vtxmasspc->GetEntry(ientry);
	vtxprob->GetEntry(ientry);
      }
      //D
      if(cat == "D"){
	vtxlen2->GetEntry(ientry);
	vtxlen12->GetEntry(ientry);
	vtxsig2->GetEntry(ientry);
	vtxsig12->GetEntry(ientry);
	vtxdirang2->GetEntry(ientry);
	vtxmult2->GetEntry(ientry);
	vtxmult->GetEntry(ientry);
	vtxmom2->GetEntry(ientry);
	vtxmass2->GetEntry(ientry);
      }
// Fill the corresponding category
      
      if(cat == "A"){
	if(nvtx_a==0){
	  trk1d0sig_histo[i]->Fill(trk1d0sig_a);
	  trk2d0sig_histo[i]->Fill(trk2d0sig_a);
	  trk1z0sig_histo[i]->Fill(trk1z0sig_a);
	  trk2z0sig_histo[i]->Fill(trk2z0sig_a);
	  trk1pt_histo[i]->Fill(trk1pt_a);
	  trk2pt_histo[i]->Fill(trk2pt_a);
	  jprobr_histo[i]->Fill(jprobr_a);
	  jprobr5sigma_histo[i]->Fill(jprobr5sigma_a);
	  jprobz_histo[i]->Fill(jprobz_a);
	  jprobz5sigma_histo[i]->Fill(jprobz5sigma_a);
	  d0bprob2_histo[i]->Fill(d0bprob2_a);
	  d0cprob2_histo[i]->Fill(d0cprob2_a);
	  d0qprob2_histo[i]->Fill(d0qprob2_a);
	  z0bprob2_histo[i]->Fill(z0bprob2_a);
	  z0cprob2_histo[i]->Fill(z0cprob2_a);
	  z0qprob2_histo[i]->Fill(z0qprob2_a);
	  nmuon_histo[i]->Fill(nmuon_a);
	  nelectron_histo[i]->Fill(nelectron_a);
	  trkmass_histo[i]->Fill(trkmass_a);
	  dEdxRatioPionOverKaonPri_histo[i]->Fill(dEdxRatioPionOverKaonPri_a);
	  dEdxRatioPionOverKaonSec_histo[i]->Fill(dEdxRatioPionOverKaonSec_a);
	}
      }
      if(cat == "B"){
        if(nvtx_a==1&&nvtxall_a==1){	  
          trk1d0sig_histo[i]->Fill(trk1d0sig_a);
          trk2d0sig_histo[i]->Fill(trk2d0sig_a);
          trk1z0sig_histo[i]->Fill(trk1z0sig_a);
          trk2z0sig_histo[i]->Fill(trk2z0sig_a);
          trk1pt_histo[i]->Fill(trk1pt_a);
          trk2pt_histo[i]->Fill(trk2pt_a);
          jprobr_histo[i]->Fill(jprobr_a);
          jprobr5sigma_histo[i]->Fill(jprobr5sigma_a);
          jprobz_histo[i]->Fill(jprobz_a);
          jprobz5sigma_histo[i]->Fill(jprobz5sigma_a);
          d0bprob2_histo[i]->Fill(d0bprob2_a);
          d0cprob2_histo[i]->Fill(d0cprob2_a);
          d0qprob2_histo[i]->Fill(d0qprob2_a);
          z0bprob2_histo[i]->Fill(z0bprob2_a);
          z0cprob2_histo[i]->Fill(z0cprob2_a);
          z0qprob2_histo[i]->Fill(z0qprob2_a);
          nmuon_histo[i]->Fill(nmuon_a);
          nelectron_histo[i]->Fill(nelectron_a);
          trkmass_histo[i]->Fill(trkmass_a);
	  dEdxRatioPionOverKaonPri_histo[i]->Fill(dEdxRatioPionOverKaonPri_a);
          dEdxRatioPionOverKaonSec_histo[i]->Fill(dEdxRatioPionOverKaonSec_a);
	  //BCD
	  one_vtxprob_histo[i]->Fill(one_vtxprob_a);
	  vtxlen1_histo[i]->Fill(vtxlen1_a);
	  vtxsig1_histo[i]->Fill(vtxsig1_a);
	  vtxdirang1_histo[i]->Fill(vtxdirang1_a);
	  vtxmult1_histo[i]->Fill(vtxmult1_a);
	  vtxmom1_histo[i]->Fill(vtxmom1_a);
	  vtxmass1_histo[i]->Fill(vtxmass1_a);
	  vtxmass_histo[i]->Fill(vtxmass_a);
	  vtxmasspc_histo[i]->Fill(vtxmasspc_a);
	  vtxprob_histo[i]->Fill(vtxprob_a);
        }
      }
      if(cat == "C"){
        if(nvtx_a==1&&nvtxall_a==2){
          trk1d0sig_histo[i]->Fill(trk1d0sig_a);
          trk2d0sig_histo[i]->Fill(trk2d0sig_a);
          trk1z0sig_histo[i]->Fill(trk1z0sig_a);
          trk2z0sig_histo[i]->Fill(trk2z0sig_a);
          trk1pt_histo[i]->Fill(trk1pt_a);
          trk2pt_histo[i]->Fill(trk2pt_a);
          jprobr_histo[i]->Fill(jprobr_a);
          jprobr5sigma_histo[i]->Fill(jprobr5sigma_a);
          jprobz_histo[i]->Fill(jprobz_a);
          jprobz5sigma_histo[i]->Fill(jprobz5sigma_a);
          d0bprob2_histo[i]->Fill(d0bprob2_a);
          d0cprob2_histo[i]->Fill(d0cprob2_a);
          d0qprob2_histo[i]->Fill(d0qprob2_a);
          z0bprob2_histo[i]->Fill(z0bprob2_a);
          z0cprob2_histo[i]->Fill(z0cprob2_a);
          z0qprob2_histo[i]->Fill(z0qprob2_a);
          nmuon_histo[i]->Fill(nmuon_a);
          nelectron_histo[i]->Fill(nelectron_a);
          trkmass_histo[i]->Fill(trkmass_a);
          dEdxRatioPionOverKaonPri_histo[i]->Fill(dEdxRatioPionOverKaonPri_a);
          dEdxRatioPionOverKaonSec_histo[i]->Fill(dEdxRatioPionOverKaonSec_a);
	  //BCD                                                                                                                                                                  
          one_vtxprob_histo[i]->Fill(one_vtxprob_a);
          vtxlen1_histo[i]->Fill(vtxlen1_a);
          vtxsig1_histo[i]->Fill(vtxsig1_a);
          vtxdirang1_histo[i]->Fill(vtxdirang1_a);
          vtxmult1_histo[i]->Fill(vtxmult1_a);
          vtxmom1_histo[i]->Fill(vtxmom1_a);
          vtxmass1_histo[i]->Fill(vtxmass1_a);
          vtxmass_histo[i]->Fill(vtxmass_a);
          vtxmasspc_histo[i]->Fill(vtxmasspc_a);
          vtxprob_histo[i]->Fill(vtxprob_a);

        }
      }
      if(cat == "D"){
        if(nvtx_a==2){
          trk1d0sig_histo[i]->Fill(trk1d0sig_a);
          trk2d0sig_histo[i]->Fill(trk2d0sig_a);
          trk1z0sig_histo[i]->Fill(trk1z0sig_a);
          trk2z0sig_histo[i]->Fill(trk2z0sig_a);
          trk1pt_histo[i]->Fill(trk1pt_a);
          trk2pt_histo[i]->Fill(trk2pt_a);
          jprobr_histo[i]->Fill(jprobr_a);
          jprobr5sigma_histo[i]->Fill(jprobr5sigma_a);
          jprobz_histo[i]->Fill(jprobz_a);
          jprobz5sigma_histo[i]->Fill(jprobz5sigma_a);
          d0bprob2_histo[i]->Fill(d0bprob2_a);
          d0cprob2_histo[i]->Fill(d0cprob2_a);
          d0qprob2_histo[i]->Fill(d0qprob2_a);
          z0bprob2_histo[i]->Fill(z0bprob2_a);
          z0cprob2_histo[i]->Fill(z0cprob2_a);
          z0qprob2_histo[i]->Fill(z0qprob2_a);
          nmuon_histo[i]->Fill(nmuon_a);
          nelectron_histo[i]->Fill(nelectron_a);
          trkmass_histo[i]->Fill(trkmass_a);
          dEdxRatioPionOverKaonPri_histo[i]->Fill(dEdxRatioPionOverKaonPri_a);
          dEdxRatioPionOverKaonSec_histo[i]->Fill(dEdxRatioPionOverKaonSec_a);
	  //BCD                                                                                                                                                                  
          one_vtxprob_histo[i]->Fill(one_vtxprob_a);
          vtxlen1_histo[i]->Fill(vtxlen1_a);
          vtxsig1_histo[i]->Fill(vtxsig1_a);
          vtxdirang1_histo[i]->Fill(vtxdirang1_a);
          vtxmult1_histo[i]->Fill(vtxmult1_a);
          vtxmom1_histo[i]->Fill(vtxmom1_a);
          vtxmass1_histo[i]->Fill(vtxmass1_a);
          vtxmass_histo[i]->Fill(vtxmass_a);
          vtxmasspc_histo[i]->Fill(vtxmasspc_a);
          vtxprob_histo[i]->Fill(vtxprob_a);
	  //D
	  vtxlen2_histo[i]->Fill(vtxlen2_a);
	  vtxlen12_histo[i]->Fill(vtxlen12_a);
	  vtxsig2_histo[i]->Fill(vtxsig2_a);
	  vtxsig12_histo[i]->Fill(vtxsig12_a);
	  vtxdirang2_histo[i]->Fill(vtxdirang2_a);
	  vtxmult2_histo[i]->Fill(vtxmult2_a);
	  vtxmult_histo[i]->Fill(vtxmult_a);
	  vtxmom2_histo[i]->Fill(vtxmom2_a);
	  vtxmass2_histo[i]->Fill(vtxmass2_a);
        }
      }
    }
    
    int entriesincat = d0bprob2_histo[i]->GetEntries();
    std::cout<<"Entries in d0bprob2 in cat. "<<cat<<": "<<entriesincat<<std::endl;
    // NORMALIZATION: 1/entriesincat
    Double_t factor=1.;
    trk1d0sig_histo[i]->Scale(factor/entriesincat);
    trk2d0sig_histo[i]->Scale(factor/entriesincat);
    trk1z0sig_histo[i]->Scale(factor/entriesincat);
    trk2z0sig_histo[i]->Scale(factor/entriesincat);
    trk1pt_histo[i]->Scale(factor/entriesincat);
    trk2pt_histo[i]->Scale(factor/entriesincat);
    jprobr_histo[i]->Scale(factor/entriesincat);
    jprobr5sigma_histo[i]->Scale(factor/entriesincat);
    jprobz_histo[i]->Scale(factor/entriesincat);
    jprobz5sigma_histo[i]->Scale(factor/entriesincat);
    d0bprob2_histo[i]->Scale(factor/entriesincat);
    d0cprob2_histo[i]->Scale(factor/entriesincat);
    d0qprob2_histo[i]->Scale(factor/entriesincat);
    z0bprob2_histo[i]->Scale(factor/entriesincat);
    z0cprob2_histo[i]->Scale(factor/entriesincat);
    z0qprob2_histo[i]->Scale(factor/entriesincat);
    nmuon_histo[i]->Scale(factor/entriesincat);
    nelectron_histo[i]->Scale(factor/entriesincat);
    trkmass_histo[i]->Scale(factor/entriesincat);
    dEdxRatioPionOverKaonPri_histo[i]->Scale(factor/entriesincat);
    dEdxRatioPionOverKaonSec_histo[i]->Scale(factor/entriesincat);
    //BCD
    if(cat == "B" || cat == "C" || cat == "D"){
      one_vtxprob_histo[i]->Scale(factor/entriesincat);
      vtxlen1_histo[i]->Scale(factor/entriesincat);
      vtxsig1_histo[i]->Scale(factor/entriesincat);
      vtxdirang1_histo[i]->Scale(factor/entriesincat);
      vtxmult1_histo[i]->Scale(factor/entriesincat);
      vtxmom1_histo[i]->Scale(factor/entriesincat);
      vtxmass1_histo[i]->Scale(factor/entriesincat);
      vtxmass_histo[i]->Scale(factor/entriesincat);
      vtxmasspc_histo[i]->Scale(factor/entriesincat);
      vtxprob_histo[i]->Scale(factor/entriesincat);
    }
    //D
    if(cat=="D"){
      vtxlen2_histo[i]->Scale(factor/entriesincat);
      vtxlen12_histo[i]->Scale(factor/entriesincat);
      vtxsig2_histo[i]->Scale(factor/entriesincat);
      vtxsig12_histo[i]->Scale(factor/entriesincat);
      vtxdirang2_histo[i]->Scale(factor/entriesincat);
      vtxmult2_histo[i]->Scale(factor/entriesincat);
      vtxmult_histo[i]->Scale(factor/entriesincat);
      vtxmom2_histo[i]->Scale(factor/entriesincat);
      vtxmass2_histo[i]->Scale(factor/entriesincat);
    }    
  }

  // PLOTS:
  
  //SetQQbarStyle();
  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  //gStyle->SetTitleBorderSize(0);
  //gStyle->SetTitleStyle(0);
  //gStyle->SetTitleX(0.2);
  //gStyle->SetMarkerSize(1.5);
  
  // d0bprob2
  d0bprob2_canvas->cd();
  //d0bprob2_canvas->SetLogy();
  d0bprob2_histo[0]->SetTitle("d0bprob2");
  d0bprob2_histo[0]->GetYaxis()->SetTitle("Entries");
  d0bprob2_histo[0]->GetXaxis()->SetTitle("d0bprob2");
  //Float_t d0bprob2_b_max = d0bprob2_histo[0]->GetMaximum();
  //Float_t d0bprob2_c_max = d0bprob2_histo[1]->GetMaximum();
  //Float_t d0bprob2_q_max = d0bprob2_histo[2]->GetMaximum();
  //Float_t d0bprob2_bc_max = TMath::Max(d0bprob2_b_max,d0bprob2_c_max);
  //Float_t d0bprob2_max = TMath::Max(d0bprob2_bc_max,d0bprob2_q_max);
  //d0bprob2_histo[0]->GetYaxis()->SetRangeUser(0.1,d0bprob2_max);
  d0bprob2_histo[0]->GetYaxis()->SetRangeUser(0,0.5);
  d0bprob2_histo[0]->Draw("histo");
  d0bprob2_histo[0]->SetLineColor(2);
  d0bprob2_histo[1]->Draw("histosame");
  d0bprob2_histo[2]->SetLineColor(3);
  d0bprob2_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
    
  TLegend *quarks_legend = new TLegend(0.6,0.2,0.8,0.3);//(0.4,0.3,0.5,0.6);                                                     
  //Fit_leg_low_W->SetTextSize(0.030);
  quarks_legend->AddEntry(d0bprob2_histo[0],"#font[42]{b quark}","l");
  quarks_legend->AddEntry(d0bprob2_histo[1],"#font[42]{c quark}","l");
  quarks_legend->AddEntry(d0bprob2_histo[2],"#font[42]{light quarks}","l");
  quarks_legend->SetFillStyle(0);                                                                                                
  quarks_legend->SetLineWidth(0);                                                                                                
  quarks_legend->SetLineColor(0);                                                                                                
  quarks_legend->SetBorderSize(0);                                                                                               
  quarks_legend->Draw();
  TString save_d0bprob2="d0bprob2_"+pol+"_cat_"+cat+".png";
  d0bprob2_canvas->SaveAs(save_d0bprob2);
  // -----------------------------------------------------------------
  // d0cprob2
  d0cprob2_canvas->cd();
  //d0cprob2_canvas->SetLogy();
  d0cprob2_histo[0]->SetTitle("d0cprob2 ");
  d0cprob2_histo[0]->GetYaxis()->SetTitle("Entries");
  d0cprob2_histo[0]->GetXaxis()->SetTitle("d0cprob2");
  /*  Float_t d0cprob2_b_max = d0cprob2_histo[0]->GetMaximum();
  Float_t d0cprob2_c_max = d0cprob2_histo[1]->GetMaximum();
  Float_t d0cprob2_q_max = d0cprob2_histo[2]->GetMaximum();
  Float_t d0cprob2_bc_max = TMath::Max(d0cprob2_b_max,d0cprob2_c_max);
  Float_t d0cprob2_max = TMath::Max(d0cprob2_bc_max,d0cprob2_q_max);
  d0cprob2_histo[0]->GetYaxis()->SetRangeUser(0.1,d0cprob2_max);
  */
  d0cprob2_histo[0]->GetYaxis()->SetRangeUser(0,0.5);
  d0cprob2_histo[0]->Draw("histo");
  d0cprob2_histo[0]->SetLineColor(2);
  d0cprob2_histo[1]->Draw("histosame");
  d0cprob2_histo[2]->SetLineColor(3);
  d0cprob2_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);  
  quarks_legend->Draw();
  TString save_d0cprob2="d0cprob2_"+pol+"_cat_"+cat+".png";
  d0cprob2_canvas->SaveAs(save_d0cprob2);
  // -----------------------------------------------------------------
  // d0qprob2
  d0qprob2_canvas->cd();
  //d0qprob2_canvas->SetLogy();
  d0qprob2_histo[0]->SetTitle("d0qprob2 ");
  d0qprob2_histo[0]->GetYaxis()->SetTitle("Entries");
  d0qprob2_histo[0]->GetXaxis()->SetTitle("d0qprob2");
  /*
  Float_t d0qprob2_b_max = d0qprob2_histo[0]->GetMaximum();
  Float_t d0qprob2_c_max = d0qprob2_histo[1]->GetMaximum();
  Float_t d0qprob2_q_max = d0qprob2_histo[2]->GetMaximum();
  Float_t d0qprob2_bc_max = TMath::Max(d0qprob2_b_max,d0qprob2_c_max);
  Float_t d0qprob2_max = TMath::Max(d0qprob2_bc_max,d0qprob2_q_max);
  d0qprob2_histo[0]->GetYaxis()->SetRangeUser(0.1,d0qprob2_max);
  */
  d0qprob2_histo[0]->GetYaxis()->SetRangeUser(0,0.5);
  d0qprob2_histo[0]->Draw("histo");
  d0qprob2_histo[0]->SetLineColor(2);
  d0qprob2_histo[1]->Draw("histosame");
  d0qprob2_histo[2]->SetLineColor(3);
  d0qprob2_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_d0qprob2="d0qprob2_"+pol+"_cat_"+cat+".png";
  d0qprob2_canvas->SaveAs(save_d0qprob2);
  // -----------------------------------------------------------------
  // z0bprob2
  z0bprob2_canvas->cd();
  //z0bprob2_canvas->SetLogy();
  z0bprob2_histo[0]->SetTitle("z0bprob2 ");
  z0bprob2_histo[0]->GetYaxis()->SetTitle("Entries");
  z0bprob2_histo[0]->GetXaxis()->SetTitle("z0bprob2");
  /*
  Float_t z0bprob2_b_max = z0bprob2_histo[0]->GetMaximum();
  Float_t z0bprob2_c_max = z0bprob2_histo[1]->GetMaximum();
  Float_t z0bprob2_q_max = z0bprob2_histo[2]->GetMaximum();
  Float_t z0bprob2_bc_max = TMath::Max(z0bprob2_b_max,z0bprob2_c_max);
  Float_t z0bprob2_max = TMath::Max(z0bprob2_bc_max,z0bprob2_q_max);
  z0bprob2_histo[0]->GetYaxis()->SetRangeUser(0.1,z0bprob2_max);
  */
  z0bprob2_histo[0]->GetYaxis()->SetRangeUser(0,0.5);
  z0bprob2_histo[0]->Draw("histo");
  z0bprob2_histo[0]->SetLineColor(2);
  z0bprob2_histo[1]->Draw("histosame");
  z0bprob2_histo[2]->SetLineColor(3);
  z0bprob2_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_z0bprob2="z0bprob2_"+pol+"_cat_"+cat+".png";
  z0bprob2_canvas->SaveAs(save_z0bprob2);
  // -----------------------------------------------------------------
  // z0cprob2
  z0cprob2_canvas->cd();
  //z0cprob2_canvas->SetLogy();
  z0cprob2_histo[0]->SetTitle("z0cprob2 ");
  z0cprob2_histo[0]->GetYaxis()->SetTitle("Entries");
  z0cprob2_histo[0]->GetXaxis()->SetTitle("z0cprob2");/*
  Float_t z0cprob2_b_max = z0cprob2_histo[0]->GetMaximum();
  Float_t z0cprob2_c_max = z0cprob2_histo[1]->GetMaximum();
  Float_t z0cprob2_q_max = z0cprob2_histo[2]->GetMaximum();
  Float_t z0cprob2_bc_max = TMath::Max(z0cprob2_b_max,z0cprob2_c_max);
  Float_t z0cprob2_max = TMath::Max(z0cprob2_bc_max,z0cprob2_q_max);
  z0cprob2_histo[0]->GetYaxis()->SetRangeUser(0.1,z0cprob2_max);*/
  z0cprob2_histo[0]->GetYaxis()->SetRangeUser(0,0.5);
  z0cprob2_histo[0]->Draw("histo");
  z0cprob2_histo[0]->SetLineColor(2);
  z0cprob2_histo[1]->Draw("histosame");
  z0cprob2_histo[2]->SetLineColor(3);
  z0cprob2_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_z0cprob2="z0cprob2_"+pol+"_cat_"+cat+".png";
  z0cprob2_canvas->SaveAs(save_z0cprob2);
  // -----------------------------------------------------------------
  // z0qprob2
  z0qprob2_canvas->cd();
  //z0qprob2_canvas->SetLogy();
  z0qprob2_histo[0]->SetTitle("z0qprob2 ");
  z0qprob2_histo[0]->GetYaxis()->SetTitle("Entries");
  z0qprob2_histo[0]->GetXaxis()->SetTitle("z0qprob2");/*
  Float_t z0qprob2_b_max = z0qprob2_histo[0]->GetMaximum();
  Float_t z0qprob2_c_max = z0qprob2_histo[1]->GetMaximum();
  Float_t z0qprob2_q_max = z0qprob2_histo[2]->GetMaximum();
  Float_t z0qprob2_bc_max = TMath::Max(z0qprob2_b_max,z0qprob2_c_max);
  Float_t z0qprob2_max = TMath::Max(z0qprob2_bc_max,z0qprob2_q_max);
  z0qprob2_histo[0]->GetYaxis()->SetRangeUser(0.1,z0qprob2_max);*/
  z0qprob2_histo[0]->GetYaxis()->SetRangeUser(0,0.5);
  z0qprob2_histo[0]->Draw("histo");
  z0qprob2_histo[0]->SetLineColor(2);
  z0qprob2_histo[1]->Draw("histosame");
  z0qprob2_histo[2]->SetLineColor(3);
  z0qprob2_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_z0qprob2="z0qprob2_"+pol+"_cat_"+cat+".png";
  z0qprob2_canvas->SaveAs(save_z0qprob2);
  // -----------------------------------------------------------------
  // trk1d0sig
  trk1d0sig_canvas->cd();
  //trk1d0sig_canvas->SetLogy();
  trk1d0sig_histo[0]->SetTitle("trk1d0sig ");
  trk1d0sig_histo[0]->GetYaxis()->SetTitle("Entries");
  trk1d0sig_histo[0]->GetXaxis()->SetTitle("trk1d0sig");/*
  Float_t trk1d0sig_b_max = trk1d0sig_histo[0]->GetMaximum();
  Float_t trk1d0sig_c_max = trk1d0sig_histo[1]->GetMaximum();
  Float_t trk1d0sig_q_max = trk1d0sig_histo[2]->GetMaximum();
  Float_t trk1d0sig_bc_max = TMath::Max(trk1d0sig_b_max,trk1d0sig_c_max);
  Float_t trk1d0sig_max = TMath::Max(trk1d0sig_bc_max,trk1d0sig_q_max);
  trk1d0sig_histo[0]->GetYaxis()->SetRangeUser(0.1,trk1d0sig_max);*/
  trk1d0sig_histo[0]->GetYaxis()->SetRangeUser(0,1);
  trk1d0sig_histo[0]->Draw("histo");
  trk1d0sig_histo[0]->SetLineColor(2);
  trk1d0sig_histo[1]->Draw("histosame");
  trk1d0sig_histo[2]->SetLineColor(3);
  trk1d0sig_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_trk1d0sig="trk1d0sig_"+pol+"_cat_"+cat+".png";
  trk1d0sig_canvas->SaveAs(save_trk1d0sig);
  // -----------------------------------------------------------------
  // trk2d0sig
  trk2d0sig_canvas->cd();
  //trk2d0sig_canvas->SetLogy();
  trk2d0sig_histo[0]->SetTitle("trk2d0sig ");
  trk2d0sig_histo[0]->GetYaxis()->SetTitle("Entries");
  trk2d0sig_histo[0]->GetXaxis()->SetTitle("trk2d0sig");/*
  Float_t trk2d0sig_b_max = trk2d0sig_histo[0]->GetMaximum();
  Float_t trk2d0sig_c_max = trk2d0sig_histo[1]->GetMaximum();
  Float_t trk2d0sig_q_max = trk2d0sig_histo[2]->GetMaximum();
  Float_t trk2d0sig_bc_max = TMath::Max(trk2d0sig_b_max,trk2d0sig_c_max);
  Float_t trk2d0sig_max = TMath::Max(trk2d0sig_bc_max,trk2d0sig_q_max);
  trk2d0sig_histo[0]->GetYaxis()->SetRangeUser(0.1,trk2d0sig_max);*/
  trk2d0sig_histo[0]->GetYaxis()->SetRangeUser(0,1);
  trk2d0sig_histo[0]->Draw("histo");
  trk2d0sig_histo[0]->SetLineColor(2);
  trk2d0sig_histo[1]->Draw("histosame");
  trk2d0sig_histo[2]->SetLineColor(3);
  trk2d0sig_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_trk2d0sig="trk2d0sig_"+pol+"_cat_"+cat+".png";
  trk2d0sig_canvas->SaveAs(save_trk2d0sig);
  // -----------------------------------------------------------------
  // trk1z0sig
  trk1z0sig_canvas->cd();
  //trk1z0sig_canvas->SetLogy();
  trk1z0sig_histo[0]->SetTitle("trk1z0sig ");
  trk1z0sig_histo[0]->GetYaxis()->SetTitle("Entries");
  trk1z0sig_histo[0]->GetXaxis()->SetTitle("trk1z0sig");
  Float_t trk1z0sig_b_max = trk1z0sig_histo[0]->GetMaximum();
  Float_t trk1z0sig_c_max = trk1z0sig_histo[1]->GetMaximum();
  Float_t trk1z0sig_q_max = trk1z0sig_histo[2]->GetMaximum();
  Float_t trk1z0sig_bc_max = TMath::Max(trk1z0sig_b_max,trk1z0sig_c_max);
  Float_t trk1z0sig_max = TMath::Max(trk1z0sig_bc_max,trk1z0sig_q_max);
  //trk1z0sig_histo[0]->GetYaxis()->SetRangeUser(0.1,trk1z0sig_max);
  trk1z0sig_histo[0]->GetYaxis()->SetRangeUser(0,trk1z0sig_max);
  trk1z0sig_histo[0]->Draw("histo");
  trk1z0sig_histo[0]->SetLineColor(2);
  trk1z0sig_histo[1]->Draw("histosame");
  trk1z0sig_histo[2]->SetLineColor(3);
  trk1z0sig_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_trk1z0sig="trk1z0sig_"+pol+"_cat_"+cat+".png";
  trk1z0sig_canvas->SaveAs(save_trk1z0sig);
  // -----------------------------------------------------------------
  // trk2z0sig
  trk2z0sig_canvas->cd();
  //trk2z0sig_canvas->SetLogy();
  trk2z0sig_histo[0]->SetTitle("trk2z0sig ");
  trk2z0sig_histo[0]->GetYaxis()->SetTitle("Entries");
  trk2z0sig_histo[0]->GetXaxis()->SetTitle("trk2z0sig");
  Float_t trk2z0sig_b_max = trk2z0sig_histo[0]->GetMaximum();
  Float_t trk2z0sig_c_max = trk2z0sig_histo[1]->GetMaximum();
  Float_t trk2z0sig_q_max = trk2z0sig_histo[2]->GetMaximum();
  Float_t trk2z0sig_bc_max = TMath::Max(trk2z0sig_b_max,trk2z0sig_c_max);
  Float_t trk2z0sig_max = TMath::Max(trk2z0sig_bc_max,trk2z0sig_q_max);
  //trk2z0sig_histo[0]->GetYaxis()->SetRangeUser(0.1,trk2z0sig_max);
  trk2z0sig_histo[0]->GetYaxis()->SetRangeUser(0,trk2z0sig_max);
  trk2z0sig_histo[0]->Draw("histo");
  trk2z0sig_histo[0]->SetLineColor(2);
  trk2z0sig_histo[1]->Draw("histosame");
  trk2z0sig_histo[2]->SetLineColor(3);
  trk2z0sig_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_trk2z0sig="trk2z0sig_"+pol+"_cat_"+cat+".png";
  trk2z0sig_canvas->SaveAs(save_trk2z0sig);
  // -----------------------------------------------------------------
  // trk1pt
  trk1pt_canvas->cd();
  //trk1pt_canvas->SetLogy();
  trk1pt_histo[0]->SetTitle("trk1pt ");
  trk1pt_histo[0]->GetYaxis()->SetTitle("Entries");
  trk1pt_histo[0]->GetXaxis()->SetTitle("trk1pt");
  Float_t trk1pt_b_max = trk1pt_histo[0]->GetMaximum();
  Float_t trk1pt_c_max = trk1pt_histo[1]->GetMaximum();
  Float_t trk1pt_q_max = trk1pt_histo[2]->GetMaximum();
  Float_t trk1pt_bc_max = TMath::Max(trk1pt_b_max,trk1pt_c_max);
  Float_t trk1pt_max = TMath::Max(trk1pt_bc_max,trk1pt_q_max);
  //trk1pt_histo[0]->GetYaxis()->SetRangeUser(0.1,trk1pt_max);
  trk1pt_histo[0]->GetYaxis()->SetRangeUser(0,trk1pt_max);
  trk1pt_histo[0]->Draw("histo");
  trk1pt_histo[0]->SetLineColor(2);
  trk1pt_histo[1]->Draw("histosame");
  trk1pt_histo[2]->SetLineColor(3);
  trk1pt_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_trk1pt="trk1pt_"+pol+"_cat_"+cat+".png";
  trk1pt_canvas->SaveAs(save_trk1pt);
  // -----------------------------------------------------------------
  // trk2pt
  trk2pt_canvas->cd();
  //trk2pt_canvas->SetLogy();
  trk2pt_histo[0]->SetTitle("trk2pt ");
  trk2pt_histo[0]->GetYaxis()->SetTitle("Entries");
  trk2pt_histo[0]->GetXaxis()->SetTitle("trk2pt");
  Float_t trk2pt_b_max = trk2pt_histo[0]->GetMaximum();
  Float_t trk2pt_c_max = trk2pt_histo[1]->GetMaximum();
  Float_t trk2pt_q_max = trk2pt_histo[2]->GetMaximum();
  Float_t trk2pt_bc_max = TMath::Max(trk2pt_b_max,trk2pt_c_max);
  Float_t trk2pt_max = TMath::Max(trk2pt_bc_max,trk2pt_q_max);
  //trk2pt_histo[0]->GetYaxis()->SetRangeUser(0.1,trk2pt_max);
  trk2pt_histo[0]->GetYaxis()->SetRangeUser(0,trk2pt_max);
  trk2pt_histo[0]->Draw("histo");
  trk2pt_histo[0]->SetLineColor(2);
  trk2pt_histo[1]->Draw("histosame");
  trk2pt_histo[2]->SetLineColor(3);
  trk2pt_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_trk2pt="trk2pt_"+pol+"_cat_"+cat+".png";
  trk2pt_canvas->SaveAs(save_trk2pt);
  // -----------------------------------------------------------------
  // jprobr
  jprobr_canvas->cd();
  //jprobr_canvas->SetLogy();
  jprobr_histo[0]->SetTitle("jprobr ");
  jprobr_histo[0]->GetYaxis()->SetTitle("Entries");
  jprobr_histo[0]->GetXaxis()->SetTitle("jprobr");
  Float_t jprobr_b_max = jprobr_histo[0]->GetMaximum();
  Float_t jprobr_c_max = jprobr_histo[1]->GetMaximum();
  Float_t jprobr_q_max = jprobr_histo[2]->GetMaximum();
  Float_t jprobr_bc_max = TMath::Max(jprobr_b_max,jprobr_c_max);
  Float_t jprobr_max = TMath::Max(jprobr_bc_max,jprobr_q_max);
  //jprobr_histo[0]->GetYaxis()->SetRangeUser(0.1,jprobr_max);
  jprobr_histo[0]->GetYaxis()->SetRangeUser(0,jprobr_max);
  jprobr_histo[0]->Draw("histo");
  jprobr_histo[0]->SetLineColor(2);
  jprobr_histo[1]->Draw("histosame");
  jprobr_histo[2]->SetLineColor(3);
  jprobr_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_jprobr="jprobr_"+pol+"_cat_"+cat+".png";
  jprobr_canvas->SaveAs(save_jprobr);
  // -----------------------------------------------------------------
  // jprobr5sigma
  jprobr5sigma_canvas->cd();
  //jprobr5sigma_canvas->SetLogy();
  jprobr5sigma_histo[0]->SetTitle("jprobr5sigma ");
  jprobr5sigma_histo[0]->GetYaxis()->SetTitle("Entries");
  jprobr5sigma_histo[0]->GetXaxis()->SetTitle("jprobr5sigma");
  Float_t jprobr5sigma_b_max = jprobr5sigma_histo[0]->GetMaximum();
  Float_t jprobr5sigma_c_max = jprobr5sigma_histo[1]->GetMaximum();
  Float_t jprobr5sigma_q_max = jprobr5sigma_histo[2]->GetMaximum();
  Float_t jprobr5sigma_bc_max = TMath::Max(jprobr5sigma_b_max,jprobr5sigma_c_max);
  Float_t jprobr5sigma_max = TMath::Max(jprobr5sigma_bc_max,jprobr5sigma_q_max);
  //jprobr5sigma_histo[0]->GetYaxis()->SetRangeUser(0.1,jprobr5sigma_max);
  jprobr5sigma_histo[0]->GetYaxis()->SetRangeUser(0,jprobr5sigma_max);
  jprobr5sigma_histo[0]->Draw("histo");
  jprobr5sigma_histo[0]->SetLineColor(2);
  jprobr5sigma_histo[1]->Draw("histosame");
  jprobr5sigma_histo[2]->SetLineColor(3);
  jprobr5sigma_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_jprobr5sigma="jprobr5sigma_"+pol+"_cat_"+cat+".png";
  jprobr5sigma_canvas->SaveAs(save_jprobr5sigma);
  // -----------------------------------------------------------------
  // jprobz
  jprobz_canvas->cd();
  //jprobz_canvas->SetLogy();
  jprobz_histo[0]->SetTitle("jprobz ");
  jprobz_histo[0]->GetYaxis()->SetTitle("Entries");
  jprobz_histo[0]->GetXaxis()->SetTitle("jprobz");
  Float_t jprobz_b_max = jprobz_histo[0]->GetMaximum();
  Float_t jprobz_c_max = jprobz_histo[1]->GetMaximum();
  Float_t jprobz_q_max = jprobz_histo[2]->GetMaximum();
  Float_t jprobz_bc_max = TMath::Max(jprobz_b_max,jprobz_c_max);
  Float_t jprobz_max = TMath::Max(jprobz_bc_max,jprobz_q_max);
  //jprobz_histo[0]->GetYaxis()->SetRangeUser(0.1,jprobz_max);
  jprobz_histo[0]->GetYaxis()->SetRangeUser(0,jprobz_max);
  jprobz_histo[0]->Draw("histo");
  jprobz_histo[0]->SetLineColor(2);
  jprobz_histo[1]->Draw("histosame");
  jprobz_histo[2]->SetLineColor(3);
  jprobz_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_jprobz="jprobz_"+pol+"_cat_"+cat+".png";
  jprobz_canvas->SaveAs(save_jprobz);
  // -----------------------------------------------------------------
  // jprobz5sigma
  jprobz5sigma_canvas->cd();
  //jprobz5sigma_canvas->SetLogy();
  jprobz5sigma_histo[0]->SetTitle("jprobz5sigma ");
  jprobz5sigma_histo[0]->GetYaxis()->SetTitle("Entries");
  jprobz5sigma_histo[0]->GetXaxis()->SetTitle("jprobz5sigma");
  Float_t jprobz5sigma_b_max = jprobz5sigma_histo[0]->GetMaximum();
  Float_t jprobz5sigma_c_max = jprobz5sigma_histo[1]->GetMaximum();
  Float_t jprobz5sigma_q_max = jprobz5sigma_histo[2]->GetMaximum();
  Float_t jprobz5sigma_bc_max = TMath::Max(jprobz5sigma_b_max,jprobz5sigma_c_max);
  Float_t jprobz5sigma_max = TMath::Max(jprobz5sigma_bc_max,jprobz5sigma_q_max);
  //jprobz5sigma_histo[0]->GetYaxis()->SetRangeUser(0.1,jprobz5sigma_max);
  jprobz5sigma_histo[0]->GetYaxis()->SetRangeUser(0,jprobz5sigma_max);
  jprobz5sigma_histo[0]->Draw("histo");
  jprobz5sigma_histo[0]->SetLineColor(2);
  jprobz5sigma_histo[1]->Draw("histosame");
  jprobz5sigma_histo[2]->SetLineColor(3);
  jprobz5sigma_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_jprobz5sigma="jprobz5sigma_"+pol+"_cat_"+cat+".png";
  jprobz5sigma_canvas->SaveAs(save_jprobz5sigma);
  // -----------------------------------------------------------------
  // nmuon
  nmuon_canvas->cd();
  //nmuon_canvas->SetLogy();
  nmuon_histo[0]->SetTitle("nmuon ");
  nmuon_histo[0]->GetYaxis()->SetTitle("Entries");
  nmuon_histo[0]->GetXaxis()->SetTitle("nmuon");
  Float_t nmuon_b_max = nmuon_histo[0]->GetMaximum();
  Float_t nmuon_c_max = nmuon_histo[1]->GetMaximum();
  Float_t nmuon_q_max = nmuon_histo[2]->GetMaximum();
  Float_t nmuon_bc_max = TMath::Max(nmuon_b_max,nmuon_c_max);
  Float_t nmuon_max = TMath::Max(nmuon_bc_max,nmuon_q_max);
  //nmuon_histo[0]->GetYaxis()->SetRangeUser(0.1,nmuon_max);
  nmuon_histo[0]->GetYaxis()->SetRangeUser(0,nmuon_max);
  nmuon_histo[0]->Draw("histo");
  nmuon_histo[0]->SetLineColor(2);
  nmuon_histo[1]->Draw("histosame");
  nmuon_histo[2]->SetLineColor(3);
  nmuon_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_nmuon="nmuon_"+pol+"_cat_"+cat+".png";
  nmuon_canvas->SaveAs(save_nmuon);
  // -----------------------------------------------------------------
  // nelectron
  nelectron_canvas->cd();
  //nelectron_canvas->SetLogy();
  nelectron_histo[0]->SetTitle("nelectron ");
  nelectron_histo[0]->GetYaxis()->SetTitle("Entries");
  nelectron_histo[0]->GetXaxis()->SetTitle("nelectron");
  Float_t nelectron_b_max = nelectron_histo[0]->GetMaximum();
  Float_t nelectron_c_max = nelectron_histo[1]->GetMaximum();
  Float_t nelectron_q_max = nelectron_histo[2]->GetMaximum();
  Float_t nelectron_bc_max = TMath::Max(nelectron_b_max,nelectron_c_max);
  Float_t nelectron_max = TMath::Max(nelectron_bc_max,nelectron_q_max);
  //nelectron_histo[0]->GetYaxis()->SetRangeUser(0.1,nelectron_max);
  nelectron_histo[0]->GetYaxis()->SetRangeUser(0,nelectron_max);
  nelectron_histo[0]->Draw("histo");
  nelectron_histo[0]->SetLineColor(2);
  nelectron_histo[1]->Draw("histosame");
  nelectron_histo[2]->SetLineColor(3);
  nelectron_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_nelectron="nelectron_"+pol+"_cat_"+cat+".png";
  nelectron_canvas->SaveAs(save_nelectron);
  // -----------------------------------------------------------------
  // trkmass
  trkmass_canvas->cd();
  //trkmass_canvas->SetLogy();
  trkmass_histo[0]->SetTitle("trkmass ");
  trkmass_histo[0]->GetYaxis()->SetTitle("Entries");
  trkmass_histo[0]->GetXaxis()->SetTitle("trkmass");
  Float_t trkmass_b_max = trkmass_histo[0]->GetMaximum();
  Float_t trkmass_c_max = trkmass_histo[1]->GetMaximum();
  Float_t trkmass_q_max = trkmass_histo[2]->GetMaximum();
  Float_t trkmass_bc_max = TMath::Max(trkmass_b_max,trkmass_c_max);
  Float_t trkmass_max = TMath::Max(trkmass_bc_max,trkmass_q_max);
  //trkmass_histo[0]->GetYaxis()->SetRangeUser(0.1,trkmass_max);
  trkmass_histo[0]->GetYaxis()->SetRangeUser(0,trkmass_max);
  trkmass_histo[0]->Draw("histo");
  trkmass_histo[0]->SetLineColor(2);
  trkmass_histo[1]->Draw("histosame");
  trkmass_histo[2]->SetLineColor(3);
  trkmass_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_trkmass="trkmass_"+pol+"_cat_"+cat+".png";
  trkmass_canvas->SaveAs(save_trkmass);
  // -----------------------
  //dEdxRatioPionOverKaonPri
  dEdxRatioPionOverKaonPri_canvas->cd();
  dEdxRatioPionOverKaonPri_histo[0]->SetTitle("dEdxRatioPionOverKaonPri ");
  dEdxRatioPionOverKaonPri_histo[0]->GetYaxis()->SetTitle("Entries (normalized)");
  dEdxRatioPionOverKaonPri_histo[0]->GetXaxis()->SetTitle("dEdxRatioPionOverKaonPri");
  Float_t dEdxRatioPionOverKaonPri_b_max = dEdxRatioPionOverKaonPri_histo[0]->GetMaximum();
  Float_t dEdxRatioPionOverKaonPri_c_max = dEdxRatioPionOverKaonPri_histo[1]->GetMaximum();
  Float_t dEdxRatioPionOverKaonPri_q_max = dEdxRatioPionOverKaonPri_histo[2]->GetMaximum();
  Float_t dEdxRatioPionOverKaonPri_bc_max = TMath::Max(dEdxRatioPionOverKaonPri_b_max,dEdxRatioPionOverKaonPri_c_max);
  Float_t dEdxRatioPionOverKaonPri_max = TMath::Max(dEdxRatioPionOverKaonPri_bc_max,dEdxRatioPionOverKaonPri_q_max);
  dEdxRatioPionOverKaonPri_histo[0]->GetYaxis()->SetRangeUser(0,1.1*dEdxRatioPionOverKaonPri_max);
  dEdxRatioPionOverKaonPri_histo[0]->Draw("histo");
  dEdxRatioPionOverKaonPri_histo[0]->SetLineColor(2);
  dEdxRatioPionOverKaonPri_histo[1]->Draw("histosame");
  dEdxRatioPionOverKaonPri_histo[2]->SetLineColor(3);
  dEdxRatioPionOverKaonPri_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_dEdxRatioPionOverKaonPri="dEdxRatioPionOverKaonPri_"+pol+"_cat_"+cat+".png";
  dEdxRatioPionOverKaonPri_canvas->SaveAs(save_dEdxRatioPionOverKaonPri);
  // ----------------------
  //dEdxRationPionOverKaonSec
  dEdxRatioPionOverKaonSec_canvas->cd();
  dEdxRatioPionOverKaonSec_histo[0]->SetTitle("dEdxRatioPionOverKaonSec ");  
  dEdxRatioPionOverKaonSec_histo[0]->GetYaxis()->SetTitle("Entries (normalized)");
  dEdxRatioPionOverKaonSec_histo[0]->GetXaxis()->SetTitle("dEdxRatioPionOverKaonSec");
  Float_t dEdxRatioPionOverKaonSec_b_max = dEdxRatioPionOverKaonSec_histo[0]->GetMaximum();
  Float_t dEdxRatioPionOverKaonSec_c_max = dEdxRatioPionOverKaonSec_histo[1]->GetMaximum();
  Float_t dEdxRatioPionOverKaonSec_q_max = dEdxRatioPionOverKaonSec_histo[2]->GetMaximum();
  Float_t dEdxRatioPionOverKaonSec_bc_max = TMath::Max(dEdxRatioPionOverKaonSec_b_max,dEdxRatioPionOverKaonSec_c_max);
  Float_t dEdxRatioPionOverKaonSec_max = TMath::Max(dEdxRatioPionOverKaonSec_bc_max,dEdxRatioPionOverKaonSec_q_max);
  dEdxRatioPionOverKaonSec_histo[0]->GetYaxis()->SetRangeUser(0,1);
  dEdxRatioPionOverKaonSec_histo[0]->Draw("histo");
  dEdxRatioPionOverKaonSec_histo[0]->SetLineColor(2);
  dEdxRatioPionOverKaonSec_histo[1]->Draw("histosame");
  dEdxRatioPionOverKaonSec_histo[2]->SetLineColor(3);
  dEdxRatioPionOverKaonSec_histo[2]->Draw("histosame");
  Labels(energy,lumi,pol,cat);
  quarks_legend->Draw();
  TString save_dEdxRatioPionOverKaonSec="dEdxRatioPionOverKaonSec_"+pol+"_cat_"+cat+".png";
  dEdxRatioPionOverKaonSec_canvas->SaveAs(save_dEdxRatioPionOverKaonSec);

  if(cat == "B" || cat == "C" || cat == "D"){
    //one_vtxprob, vtxlen1, vtxsig1, vtxmult1, vtxmom1, vtxmass1, vtxmass, vtxmasspc, vtxprob
    // -----------------------------------------------------------------
    // one_vtxprob
    one_vtxprob_canvas->cd();
    //one_vtxprob_canvas->SetLogy();
    one_vtxprob_histo[0]->SetTitle("one_vtxprob ");
    one_vtxprob_histo[0]->GetYaxis()->SetTitle("Entries");
    one_vtxprob_histo[0]->GetXaxis()->SetTitle("one_vtxprob");
    Float_t one_vtxprob_b_max = one_vtxprob_histo[0]->GetMaximum();
    Float_t one_vtxprob_c_max = one_vtxprob_histo[1]->GetMaximum();
    Float_t one_vtxprob_q_max = one_vtxprob_histo[2]->GetMaximum();
    Float_t one_vtxprob_bc_max = TMath::Max(one_vtxprob_b_max,one_vtxprob_c_max);
    Float_t one_vtxprob_max = TMath::Max(one_vtxprob_bc_max,one_vtxprob_q_max);
    //one_vtxprob_histo[0]->GetYaxis()->SetRangeUser(0.1,one_vtxprob_max);
    one_vtxprob_histo[0]->GetYaxis()->SetRangeUser(0,one_vtxprob_max);
    one_vtxprob_histo[0]->Draw("histo");
    one_vtxprob_histo[0]->SetLineColor(2);
    one_vtxprob_histo[1]->Draw("histosame");
    one_vtxprob_histo[2]->SetLineColor(3);
    one_vtxprob_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_one_vtxprob="one_vtxprob_"+pol+"_cat_"+cat+".png";
    one_vtxprob_canvas->SaveAs(save_one_vtxprob);
    // -----------------------------------------------------------------
    // vtxlen1
    vtxlen1_canvas->cd();
    //vtxlen1_canvas->SetLogy();
    vtxlen1_histo[0]->SetTitle("vtxlen1 ");
    vtxlen1_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxlen1_histo[0]->GetXaxis()->SetTitle("vtxlen1");
    Float_t vtxlen1_b_max = vtxlen1_histo[0]->GetMaximum();
    Float_t vtxlen1_c_max = vtxlen1_histo[1]->GetMaximum();
    Float_t vtxlen1_q_max = vtxlen1_histo[2]->GetMaximum();
    Float_t vtxlen1_bc_max = TMath::Max(vtxlen1_b_max,vtxlen1_c_max);
    Float_t vtxlen1_max = TMath::Max(vtxlen1_bc_max,vtxlen1_q_max);
    //vtxlen1_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxlen1_max);
    vtxlen1_histo[0]->GetYaxis()->SetRangeUser(0,vtxlen1_max);
    vtxlen1_histo[0]->Draw("histo");
    vtxlen1_histo[0]->SetLineColor(2);
    vtxlen1_histo[1]->Draw("histosame");
    vtxlen1_histo[2]->SetLineColor(3);
    vtxlen1_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxlen1="vtxlen1_"+pol+"_cat_"+cat+".png";
    vtxlen1_canvas->SaveAs(save_vtxlen1);
    // -----------------------------------------------------------------
    // vtxsig1
    vtxsig1_canvas->cd();
    //vtxsig1_canvas->SetLogy();
    vtxsig1_histo[0]->SetTitle("vtxsig1 ");
    vtxsig1_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxsig1_histo[0]->GetXaxis()->SetTitle("vtxsig1");
    Float_t vtxsig1_b_max = vtxsig1_histo[0]->GetMaximum();
    Float_t vtxsig1_c_max = vtxsig1_histo[1]->GetMaximum();
    Float_t vtxsig1_q_max = vtxsig1_histo[2]->GetMaximum();
    Float_t vtxsig1_bc_max = TMath::Max(vtxsig1_b_max,vtxsig1_c_max);
    Float_t vtxsig1_max = TMath::Max(vtxsig1_bc_max,vtxsig1_q_max);
    //vtxsig1_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxsig1_max);
    vtxsig1_histo[0]->GetYaxis()->SetRangeUser(0,vtxsig1_max);
    vtxsig1_histo[0]->Draw("histo");
    vtxsig1_histo[0]->SetLineColor(2);
    vtxsig1_histo[1]->Draw("histosame");
    vtxsig1_histo[2]->SetLineColor(3);
    vtxsig1_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxsig1="vtxsig1_"+pol+"_cat_"+cat+".png";
    vtxsig1_canvas->SaveAs(save_vtxsig1);
    // -----------------------------------------------------------------
    // vtxdirang1
    vtxdirang1_canvas->cd();
    //vtxdirang1_canvas->SetLogy();
    vtxdirang1_histo[0]->SetTitle("vtxdirang1 ");
    vtxdirang1_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxdirang1_histo[0]->GetXaxis()->SetTitle("vtxdirang1");
    Float_t vtxdirang1_b_max = vtxdirang1_histo[0]->GetMaximum();
    Float_t vtxdirang1_c_max = vtxdirang1_histo[1]->GetMaximum();
    Float_t vtxdirang1_q_max = vtxdirang1_histo[2]->GetMaximum();
    Float_t vtxdirang1_bc_max = TMath::Max(vtxdirang1_b_max,vtxdirang1_c_max);
    Float_t vtxdirang1_max = TMath::Max(vtxdirang1_bc_max,vtxdirang1_q_max);
    //vtxdirang1_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxdirang1_max);
    vtxdirang1_histo[0]->GetYaxis()->SetRangeUser(0,vtxdirang1_max);
    vtxdirang1_histo[0]->Draw("histo");
    vtxdirang1_histo[0]->SetLineColor(2);
    vtxdirang1_histo[1]->Draw("histosame");
    vtxdirang1_histo[2]->SetLineColor(3);
    vtxdirang1_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxdirang1="vtxdirang1_"+pol+"_cat_"+cat+".png";
    vtxdirang1_canvas->SaveAs(save_vtxdirang1);
    // -----------------------------------------------------------------
    // vtxmult1
    vtxmult1_canvas->cd();
    //vtxmult1_canvas->SetLogy();
    vtxmult1_histo[0]->SetTitle("vtxmult1 ");
    vtxmult1_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxmult1_histo[0]->GetXaxis()->SetTitle("vtxmult1");
    Float_t vtxmult1_b_max = vtxmult1_histo[0]->GetMaximum();
    Float_t vtxmult1_c_max = vtxmult1_histo[1]->GetMaximum();
    Float_t vtxmult1_q_max = vtxmult1_histo[2]->GetMaximum();
    Float_t vtxmult1_bc_max = TMath::Max(vtxmult1_b_max,vtxmult1_c_max);
    Float_t vtxmult1_max = TMath::Max(vtxmult1_bc_max,vtxmult1_q_max);
    //vtxmult1_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxmult1_max);
    vtxmult1_histo[0]->GetYaxis()->SetRangeUser(0,vtxmult1_max);
    vtxmult1_histo[0]->Draw("histo");
    vtxmult1_histo[0]->SetLineColor(2);
    vtxmult1_histo[1]->Draw("histosame");
    vtxmult1_histo[2]->SetLineColor(3);
    vtxmult1_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxmult1="vtxmult1_"+pol+"_cat_"+cat+".png";
    vtxmult1_canvas->SaveAs(save_vtxmult1);
    // -----------------------------------------------------------------
    // vtxmom1
    vtxmom1_canvas->cd();
    //vtxmom1_canvas->SetLogy();
    vtxmom1_histo[0]->SetTitle("vtxmom1 ");
    vtxmom1_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxmom1_histo[0]->GetXaxis()->SetTitle("vtxmom1");
    Float_t vtxmom1_b_max = vtxmom1_histo[0]->GetMaximum();
    Float_t vtxmom1_c_max = vtxmom1_histo[1]->GetMaximum();
    Float_t vtxmom1_q_max = vtxmom1_histo[2]->GetMaximum();
    Float_t vtxmom1_bc_max = TMath::Max(vtxmom1_b_max,vtxmom1_c_max);
    Float_t vtxmom1_max = TMath::Max(vtxmom1_bc_max,vtxmom1_q_max);
    //vtxmom1_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxmom1_max);
    vtxmom1_histo[0]->GetYaxis()->SetRangeUser(0,vtxmom1_max);
    vtxmom1_histo[0]->Draw("histo");
    vtxmom1_histo[0]->SetLineColor(2);
    vtxmom1_histo[1]->Draw("histosame");
    vtxmom1_histo[2]->SetLineColor(3);
    vtxmom1_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxmom1="vtxmom1_"+pol+"_cat_"+cat+".png";
    vtxmom1_canvas->SaveAs(save_vtxmom1);
    // -----------------------------------------------------------------
    // vtxmass1
    vtxmass1_canvas->cd();
    //vtxmass1_canvas->SetLogy();
    vtxmass1_histo[0]->SetTitle("vtxmass1 ");
    vtxmass1_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxmass1_histo[0]->GetXaxis()->SetTitle("vtxmass");
    Float_t vtxmass1_b_max = vtxmass1_histo[0]->GetMaximum();
    Float_t vtxmass1_c_max = vtxmass1_histo[1]->GetMaximum();
    Float_t vtxmass1_q_max = vtxmass1_histo[2]->GetMaximum();
    Float_t vtxmass1_bc_max = TMath::Max(vtxmass1_b_max,vtxmass1_c_max);
    Float_t vtxmass1_max = TMath::Max(vtxmass1_bc_max,vtxmass1_q_max);
    //vtxmass1_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxmass1_max);
    vtxmass1_histo[0]->GetYaxis()->SetRangeUser(0,vtxmass1_max);
    vtxmass1_histo[0]->Draw("histo");
    vtxmass1_histo[0]->SetLineColor(2);
    vtxmass1_histo[1]->Draw("histosame");
    vtxmass1_histo[2]->SetLineColor(3);
    vtxmass1_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxmass1="vtxmass1_"+pol+"_cat_"+cat+".png";
    vtxmass1_canvas->SaveAs(save_vtxmass1);
    // -----------------------------------------------------------------
    // vtxmass                                                                                                                   
    vtxmass_canvas->cd();
    //vtxmass_canvas->SetLogy();
    vtxmass_histo[0]->SetTitle("vtxmass ");
    vtxmass_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxmass_histo[0]->GetXaxis()->SetTitle("vtxmass");
    Float_t vtxmass_b_max = vtxmass1_histo[0]->GetMaximum();
    Float_t vtxmass_c_max = vtxmass1_histo[1]->GetMaximum();
    Float_t vtxmass_q_max = vtxmass1_histo[2]->GetMaximum();
    Float_t vtxmass_bc_max = TMath::Max(vtxmass_b_max,vtxmass_c_max);
    Float_t vtxmass_max = TMath::Max(vtxmass_bc_max,vtxmass_q_max);
    //vtxmass_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxmass_max);
    vtxmass_histo[0]->GetYaxis()->SetRangeUser(0,vtxmass_max);
    vtxmass_histo[0]->Draw("histo");
    vtxmass_histo[0]->SetLineColor(2);
    vtxmass_histo[1]->Draw("histosame");
    vtxmass_histo[2]->SetLineColor(3);
    vtxmass_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxmass="vtxmass_"+pol+"_cat_"+cat+".png";
    vtxmass_canvas->SaveAs(save_vtxmass);
    // -----------------------------------------------------------------
    // vtxmasspc
    vtxmasspc_canvas->cd();
    //vtxmasspc_canvas->SetLogy();
    vtxmasspc_histo[0]->SetTitle("vtxmasspc ");
    vtxmasspc_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxmasspc_histo[0]->GetXaxis()->SetTitle("vtxmasspc");
    Float_t vtxmasspc_b_max = vtxmasspc_histo[0]->GetMaximum();
    Float_t vtxmasspc_c_max = vtxmasspc_histo[1]->GetMaximum();
    Float_t vtxmasspc_q_max = vtxmasspc_histo[2]->GetMaximum();
    Float_t vtxmasspc_bc_max = TMath::Max(vtxmasspc_b_max,vtxmasspc_c_max);
    Float_t vtxmasspc_max = TMath::Max(vtxmasspc_bc_max,vtxmasspc_q_max);
    //vtxmasspc_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxmasspc_max);
    vtxmasspc_histo[0]->GetYaxis()->SetRangeUser(0,vtxmasspc_max);
    vtxmasspc_histo[0]->Draw("histo");
    vtxmasspc_histo[0]->SetLineColor(2);
    vtxmasspc_histo[1]->Draw("histosame");
    vtxmasspc_histo[2]->SetLineColor(3);
    vtxmasspc_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxmasspc="vtxmasspc_"+pol+"_cat_"+cat+".png";
    vtxmasspc_canvas->SaveAs(save_vtxmasspc);
    // -----------------------------------------------------------------
    // vtxprob
    vtxprob_canvas->cd();
    //vtxprob_canvas->SetLogy();
    vtxprob_histo[0]->SetTitle("vtxprob ");
    vtxprob_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxprob_histo[0]->GetXaxis()->SetTitle("vtxprob");
    Float_t vtxprob_b_max = vtxprob_histo[0]->GetMaximum();
    Float_t vtxprob_c_max = vtxprob_histo[1]->GetMaximum();
    Float_t vtxprob_q_max = vtxprob_histo[2]->GetMaximum();
    Float_t vtxprob_bc_max = TMath::Max(vtxprob_b_max,vtxprob_c_max);
    Float_t vtxprob_max = TMath::Max(vtxprob_bc_max,vtxprob_q_max);
    //vtxprob_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxprob_max);
    vtxprob_histo[0]->GetYaxis()->SetRangeUser(0,vtxprob_max);
    vtxprob_histo[0]->Draw("histo");
    vtxprob_histo[0]->SetLineColor(2);
    vtxprob_histo[1]->Draw("histosame");
    vtxprob_histo[2]->SetLineColor(3);
    vtxprob_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxprob="vtxprob_"+pol+"_cat_"+cat+".png";
    vtxprob_canvas->SaveAs(save_vtxprob);
  }

  if(cat == "D"){
    //vtxlen2, vtxlen12, vtxsig2, vtxsig12, vtxdirang2, vtxmult2, vtxmult, vtxmom2, vtxmass2
    // -----------------------------------------------------------------
    // vtxlen2
    vtxlen2_canvas->cd();
    //vtxlen2_canvas->SetLogy();
    vtxlen2_histo[0]->SetTitle("vtxlen2 ");
    vtxlen2_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxlen2_histo[0]->GetXaxis()->SetTitle("vtxlen2");
    Float_t vtxlen2_b_max = vtxlen2_histo[0]->GetMaximum();
    Float_t vtxlen2_c_max = vtxlen2_histo[1]->GetMaximum();
    Float_t vtxlen2_q_max = vtxlen2_histo[2]->GetMaximum();
    Float_t vtxlen2_bc_max = TMath::Max(vtxlen2_b_max,vtxlen2_c_max);
    Float_t vtxlen2_max = TMath::Max(vtxlen2_bc_max,vtxlen2_q_max);
    //vtxlen2_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxlen2_max);
    vtxlen2_histo[0]->GetYaxis()->SetRangeUser(0,vtxlen2_max);
    vtxlen2_histo[0]->Draw("histo");
    vtxlen2_histo[0]->SetLineColor(2);
    vtxlen2_histo[1]->Draw("histosame");
    vtxlen2_histo[2]->SetLineColor(3);
    vtxlen2_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxlen2="vtxlen2_"+pol+"_cat_"+cat+".png";
    vtxlen2_canvas->SaveAs(save_vtxlen2);
    // -----------------------------------------------------------------                                                                                      
    // vtxlen12
    vtxlen12_canvas->cd();
    //vtxlen12_canvas->SetLogy();
    vtxlen12_histo[0]->SetTitle("vtxlen12 ");
    vtxlen12_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxlen12_histo[0]->GetXaxis()->SetTitle("vtxlen12");
    Float_t vtxlen12_b_max = vtxlen12_histo[0]->GetMaximum();
    Float_t vtxlen12_c_max = vtxlen12_histo[1]->GetMaximum();
    Float_t vtxlen12_q_max = vtxlen12_histo[2]->GetMaximum();
    Float_t vtxlen12_bc_max = TMath::Max(vtxlen12_b_max,vtxlen12_c_max);
    Float_t vtxlen12_max = TMath::Max(vtxlen12_bc_max,vtxlen12_q_max);
    //vtxlen12_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxlen12_max);
    vtxlen12_histo[0]->GetYaxis()->SetRangeUser(0,vtxlen12_max);
    vtxlen12_histo[0]->Draw("histo");
    vtxlen12_histo[0]->SetLineColor(2);
    vtxlen12_histo[1]->Draw("histosame");
    vtxlen12_histo[2]->SetLineColor(3);
    vtxlen12_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxlen12="vtxlen12_"+pol+"_cat_"+cat+".png";
    vtxlen12_canvas->SaveAs(save_vtxlen12);
    // -----------------------------------------------------------------
    // vtxsig2
    vtxsig2_canvas->cd();
    //vtxsig2_canvas->SetLogy();
    vtxsig2_histo[0]->SetTitle("vtxsig2 ");
    vtxsig2_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxsig2_histo[0]->GetXaxis()->SetTitle("vtxsig2");
    Float_t vtxsig2_b_max = vtxsig2_histo[0]->GetMaximum();
    Float_t vtxsig2_c_max = vtxsig2_histo[1]->GetMaximum();
    Float_t vtxsig2_q_max = vtxsig2_histo[2]->GetMaximum();
    Float_t vtxsig2_bc_max = TMath::Max(vtxsig2_b_max,vtxsig2_c_max);
    Float_t vtxsig2_max = TMath::Max(vtxsig2_bc_max,vtxsig2_q_max);
    //vtxsig2_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxsig2_max);
    vtxsig2_histo[0]->GetYaxis()->SetRangeUser(0,vtxsig2_max);
    vtxsig2_histo[0]->Draw("histo");
    vtxsig2_histo[0]->SetLineColor(2);
    vtxsig2_histo[1]->Draw("histosame");
    vtxsig2_histo[2]->SetLineColor(3);
    vtxsig2_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxsig2="vtxsig2_"+pol+"_cat_"+cat+".png";
    vtxsig2_canvas->SaveAs(save_vtxsig2);
    // -----------------------------------------------------------------
    // vtxsig12
    vtxsig12_canvas->cd();
    //vtxsig12_canvas->SetLogy();
    vtxsig12_histo[0]->SetTitle("vtxsig12 ");
    vtxsig12_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxsig12_histo[0]->GetXaxis()->SetTitle("vtxsig12");
    Float_t vtxsig12_b_max = vtxsig12_histo[0]->GetMaximum();
    Float_t vtxsig12_c_max = vtxsig12_histo[1]->GetMaximum();
    Float_t vtxsig12_q_max = vtxsig12_histo[2]->GetMaximum();
    Float_t vtxsig12_bc_max = TMath::Max(vtxsig12_b_max,vtxsig12_c_max);
    Float_t vtxsig12_max = TMath::Max(vtxsig12_bc_max,vtxsig12_q_max);
    //vtxsig12_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxsig12_max);
    vtxsig12_histo[0]->GetYaxis()->SetRangeUser(0,vtxsig12_max);
    vtxsig12_histo[0]->Draw("histo");
    vtxsig12_histo[0]->SetLineColor(2);
    vtxsig12_histo[1]->Draw("histosame");
    vtxsig12_histo[2]->SetLineColor(3);
    vtxsig12_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxsig12="vtxsig12_"+pol+"_cat_"+cat+".png";
    vtxsig12_canvas->SaveAs(save_vtxsig12);
    // -----------------------------------------------------------------
    // vtxdirang2
    vtxdirang2_canvas->cd();
    //vtxdirang2_canvas->SetLogy();
    vtxdirang2_histo[0]->SetTitle("vtxdirang2 ");
    vtxdirang2_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxdirang2_histo[0]->GetXaxis()->SetTitle("vtxdirang2");
    Float_t vtxdirang2_b_max = vtxdirang2_histo[0]->GetMaximum();
    Float_t vtxdirang2_c_max = vtxdirang2_histo[1]->GetMaximum();
    Float_t vtxdirang2_q_max = vtxdirang2_histo[2]->GetMaximum();
    Float_t vtxdirang2_bc_max = TMath::Max(vtxdirang2_b_max,vtxdirang2_c_max);
    Float_t vtxdirang2_max = TMath::Max(vtxdirang2_bc_max,vtxdirang2_q_max);
    //vtxdirang2_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxdirang2_max);
    vtxdirang2_histo[0]->GetYaxis()->SetRangeUser(0,vtxdirang2_max);
    vtxdirang2_histo[0]->Draw("histo");
    vtxdirang2_histo[0]->SetLineColor(2);
    vtxdirang2_histo[1]->Draw("histosame");
    vtxdirang2_histo[2]->SetLineColor(3);
    vtxdirang2_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxdirang2="vtxdirang2_"+pol+"_cat_"+cat+".png";
    vtxdirang2_canvas->SaveAs(save_vtxdirang2);
    // -----------------------------------------------------------------
    // vtxmult2
    vtxmult2_canvas->cd();
    //vtxmult2_canvas->SetLogy();
    vtxmult2_histo[0]->SetTitle("vtxmult2 ");
    vtxmult2_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxmult2_histo[0]->GetXaxis()->SetTitle("vtxmult2");
    Float_t vtxmult2_b_max = vtxmult2_histo[0]->GetMaximum();
    Float_t vtxmult2_c_max = vtxmult2_histo[1]->GetMaximum();
    Float_t vtxmult2_q_max = vtxmult2_histo[2]->GetMaximum();
    Float_t vtxmult2_bc_max = TMath::Max(vtxmult2_b_max,vtxmult2_c_max);
    Float_t vtxmult2_max = TMath::Max(vtxmult2_bc_max,vtxmult2_q_max);
    //vtxmult2_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxmult2_max);
    vtxmult2_histo[0]->GetYaxis()->SetRangeUser(0,vtxmult2_max);
    vtxmult2_histo[0]->Draw("histo");
    vtxmult2_histo[0]->SetLineColor(2);
    vtxmult2_histo[1]->Draw("histosame");
    vtxmult2_histo[2]->SetLineColor(3);
    vtxmult2_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxmult2="vtxmult2_"+pol+"_cat_"+cat+".png";
    vtxmult2_canvas->SaveAs(save_vtxmult2);
    // -----------------------------------------------------------------
    // vtxmult
    vtxmult_canvas->cd();
    //vtxmult_canvas->SetLogy();
    vtxmult_histo[0]->SetTitle("vtxmult ");
    vtxmult_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxmult_histo[0]->GetXaxis()->SetTitle("vtxmult");
    Float_t vtxmult_b_max = vtxmult_histo[0]->GetMaximum();
    Float_t vtxmult_c_max = vtxmult_histo[1]->GetMaximum();
    Float_t vtxmult_q_max = vtxmult_histo[2]->GetMaximum();
    Float_t vtxmult_bc_max = TMath::Max(vtxmult_b_max,vtxmult_c_max);
    Float_t vtxmult_max = TMath::Max(vtxmult_bc_max,vtxmult_q_max);
    //vtxmult_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxmult_max);
    vtxmult_histo[0]->GetYaxis()->SetRangeUser(0,vtxmult_max);
    vtxmult_histo[0]->Draw("histo");
    vtxmult_histo[0]->SetLineColor(2);
    vtxmult_histo[1]->Draw("histosame");
    vtxmult_histo[2]->SetLineColor(3);
    vtxmult_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxmult="vtxmult_"+pol+"_cat_"+cat+".png";
    vtxmult_canvas->SaveAs(save_vtxmult);
    // -----------------------------------------------------------------
    // vtxmom2
    vtxmom2_canvas->cd();
    //vtxmom2_canvas->SetLogy();
    vtxmom2_histo[0]->SetTitle("vtxmom2 ");
    vtxmom2_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxmom2_histo[0]->GetXaxis()->SetTitle("vtxmom2");
    Float_t vtxmom2_b_max = vtxmom2_histo[0]->GetMaximum();
    Float_t vtxmom2_c_max = vtxmom2_histo[1]->GetMaximum();
    Float_t vtxmom2_q_max = vtxmom2_histo[2]->GetMaximum();
    Float_t vtxmom2_bc_max = TMath::Max(vtxmom2_b_max,vtxmom2_c_max);
    Float_t vtxmom2_max = TMath::Max(vtxmom2_bc_max,vtxmom2_q_max);
    //vtxmom2_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxmom2_max);
    vtxmom2_histo[0]->GetYaxis()->SetRangeUser(0,vtxmom2_max);
    vtxmom2_histo[0]->Draw("histo");
    vtxmom2_histo[0]->SetLineColor(2);
    vtxmom2_histo[1]->Draw("histosame");
    vtxmom2_histo[2]->SetLineColor(3);
    vtxmom2_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxmom2="vtxmom2_"+pol+"_cat_"+cat+".png";
    vtxmom2_canvas->SaveAs(save_vtxmom2);
    // -----------------------------------------------------------------
    // vtxmass2
    vtxmass2_canvas->cd();
    //vtxmass2_canvas->SetLogy();
    vtxmass2_histo[0]->SetTitle("vtxmass2 ");
    vtxmass2_histo[0]->GetYaxis()->SetTitle("Entries");
    vtxmass2_histo[0]->GetXaxis()->SetTitle("vtxmass2");
    Float_t vtxmass2_b_max = vtxmass2_histo[0]->GetMaximum();
    Float_t vtxmass2_c_max = vtxmass2_histo[1]->GetMaximum();
    Float_t vtxmass2_q_max = vtxmass2_histo[2]->GetMaximum();
    Float_t vtxmass2_bc_max = TMath::Max(vtxmass2_b_max,vtxmass2_c_max);
    Float_t vtxmass2_max = TMath::Max(vtxmass2_bc_max,vtxmass2_q_max);
    //vtxmass2_histo[0]->GetYaxis()->SetRangeUser(0.1,vtxmass2_max);
    vtxmass2_histo[0]->GetYaxis()->SetRangeUser(0,vtxmass2_max);
    vtxmass2_histo[0]->Draw("histo");
    vtxmass2_histo[0]->SetLineColor(2);
    vtxmass2_histo[1]->Draw("histosame");
    vtxmass2_histo[2]->SetLineColor(3);
    vtxmass2_histo[2]->Draw("histosame");
    Labels(energy,lumi,pol,cat);
    quarks_legend->Draw();
    TString save_vtxmass2="vtxmass2_"+pol+"_cat_"+cat+".png";
    vtxmass2_canvas->SaveAs(save_vtxmass2);
  }
  
}
