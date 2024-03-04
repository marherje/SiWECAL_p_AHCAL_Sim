#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <numeric>

using namespace std;

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

#define N_ENERGIES 18 //18

void get_res(int &nhit,     int &nhit_masked,
             float &sume,   float &sume_masked,
             float &weight, float &weight_masked,
             vector<float> * hit_energy, vector<int> *hit_slab, vector<int> *hit_isMasked,
            //  double (&W_thicknesses)[15]) {
             TVectorD W_thicknesses) {
  if(hit_energy->size() > 0){    
    //cout<<W_thicknesses.Min()<<endl;
    // First option: use the minimum of the used
    // Second option: use the paper as reference 0.4X0, X0=3.5mm
   for (int j = 0; j < hit_energy->size(); j++) {
            if (hit_isMasked->at(j) == 0) {
                nhit_masked += 1;
                sume_masked += hit_energy->at(j);
                weight_masked += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
            }
            nhit += 1;
            sume += hit_energy->at(j);
            weight += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
        }
  }
    return;
}

double hit_layer(int layer, vector<int> *hit_slab, vector<int> *hit_isMasked, bool masked=false, bool normalized=false){
  double count = 0.;
  double weight = 1.;

  if(hit_slab->size() > 0){
    for( int j = 0; j < hit_slab->size(); j++ ) {
      if( masked && hit_isMasked->at(j) == 1 ) continue;
      if( hit_slab->at(j) == layer ) count += 1;
    }
    if(normalized == true) weight = 1./hit_slab->size();
  }
  return weight*count;
}

float moliere(vector<float> * hit_energy,
              vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z,
              vector<int> * hit_isMasked, bool masked=false, float containment = 0.90) {

    vector<float> hit_rs;
    vector<float> hit_es;
    float mol_rad = 0.;
    float sume = 0.;
    float wx = 0.; float wy = 0.; float wz = 0.;
    float r;
    if(hit_energy->size() > 0){
      //cout<<"before first loop in molfunc"<<endl;
    
      for (int j = 0; j < hit_energy->size(); j++){
        if (masked && hit_isMasked->at(j) == 1) continue;
        sume += hit_energy->at(j);
        wx += hit_x->at(j) * hit_energy->at(j);
        wy += hit_y->at(j) * hit_energy->at(j);
        wz += hit_z->at(j) * hit_energy->at(j);
      }
      
      float bary_x = wx / sume;  
      float bary_y = wy / sume;  
      if( (sume == 0) or (sume < 0) ){
	bary_x = 0;
	bary_y = 0;
      }
      //cout<<"bary_x: "<<bary_x<<", bary_y: "<<bary_y<<endl;
      
      //cout<<"before second loop in molfunc"<<endl;
      for (int j = 0; j < hit_energy->size(); j++){
        if (masked && hit_isMasked->at(j) == 1) continue;
        r = pow(pow((hit_x->at(j) - bary_x) , 2) + pow((hit_y->at(j) - bary_y), 2), 0.5);
        hit_rs.push_back(r);
      }
      for (auto k: sort_indexes(hit_rs)) hit_es.push_back(hit_energy->at(k));
      
      sort(hit_rs.begin(), hit_rs.end());
      
      float mol_e = 0.;
      int mol_i=0;    
      //cout<<"before third loop in molfunc"<<endl;
      for (int j = 0; j < hit_rs.size(); j++) {
	mol_e += hit_es.at(j);
	if (mol_e >= containment * sume){
	  mol_i=j-1;
	  break; 
	}
      }
      
      if(mol_i<0) mol_i=0;
      mol_rad=hit_rs.at(mol_i);
    }
    return mol_rad;      

}


void fit_res(TH1 * h, double &mu, double &sig, double &res){
    h->Fit("gaus", "q");
    mu = h->GetFunction("gaus")->GetParameter("Mean");
    sig = h->GetFunction("gaus")->GetParameter("Sigma");
    res = h->GetFunction("gaus")->GetParameter("Sigma") / h->GetFunction("gaus")->GetParameter("Mean");
    return;
}

void fit_mol(TH1 * h, double &G_mu, double &G_sig, double &L_mpv, double &L_fwhm, double &MPV, double &FWHM, double &xmin, double &xmax){

  /*
  TF1 *GaussDist = new TF1("GaussDist","[0]*exp((-0.5*((x-[1])/[2])^2))/([2]*sqrt(2*3.14159265358979323846))",5.,100.);
  GaussDist->SetParameter(0,mu*5000);
  GaussDist->SetParameter(1,mu);
  GaussDist->SetParameter(2,sig);
  
  TF1 *LandauDist = new TF1("LandauDist","TMath::Landau(x,[3],[4])",5.,100.);
  LandauDist->SetParameter(3,mpv);
  LandauDist->SetParameter(4,fwhm);
  */

  int h_binmax = h->GetMaximumBin();
  int h_NBINS = h->GetNbinsX();
  //double h_xmaxvalue = h->GetXaxis()->GetBinCenter(h_binmax);
  double h_MEAN = h->GetMean();
  double h_maxvalue = h->GetBinContent(h_binmax);
  double h_RMS = h->GetRMS();

  TF1 *LandauGaussDist = new TF1("LandauGaussDist","[0]*exp((-0.5*((x-[1])/[2])^2))/([2]*sqrt(2*TMath::Pi()))*TMath::Landau(x,[3],[4])",xmin,xmax);
  cout<<"Par. mol.: "<<0.5*h_NBINS*h_maxvalue<<", "<<h_MEAN<<", "<<h_RMS<<", "<<h_MEAN<<", "<<h_RMS<<endl;
  LandauGaussDist->SetParameter(0,0.5*h_NBINS*h_maxvalue);
  LandauGaussDist->SetParameter(1,0.);
  LandauGaussDist->SetParameter(2,h_RMS);
  LandauGaussDist->SetParameter(3,h_MEAN);
  LandauGaussDist->SetParameter(4,h_RMS);

  h->Fit(LandauGaussDist,"LRq");
  G_mu = LandauGaussDist->GetParameter(1);
  G_sig = LandauGaussDist->GetParameter(2);
  L_mpv = LandauGaussDist->GetParameter(3);
  L_fwhm = LandauGaussDist->GetParameter(4);

  MPV = LandauGaussDist->GetMaximumX();
  
  double maxvalue = LandauGaussDist->Eval(MPV);
  cout<<"Max. Value: "<<maxvalue<<", at x= "<<MPV<<endl;
  double x1=0.;
  double x2=MPV;
  for(int i=0;i<10000000;i++){
    x1 += 0.00001;
    double value = LandauGaussDist->Eval(x1);
    double distance = abs(maxvalue/2 - value);
    //cout<<distance<<endl;
    if(distance < 0.01) break;
  }
  for(int i=0;i<10000000;i++){
    x2 += 0.00001;
    double value = LandauGaussDist->Eval(x2);
    double distance = abs(maxvalue/2 - value);
    //cout<<distance<<endl;
    if(abs(maxvalue/2-value) < 0.01) break;
  }
  
  cout<<"x1: "<<x1<<", x2: "<<x2<<endl;
  FWHM = x2 - x1;
  
  return; 
}

void fit_shower(TH1 * h, double &G_mu, double &G_sig, double &L_mpv, double &L_fwhm, double &MPV, double &FWHM, double &xmin, double &xmax){

  int h_binmax = h->GetMaximumBin();
  int h_NBINS = h->GetNbinsX();
  //double h_xmaxvalue = h->GetXaxis()->GetBinCenter(h_binmax);
  double h_MEAN= h->GetMean();
  double h_maxvalue = h->GetBinContent(h_binmax);
  double h_RMS = h->GetRMS();
  
  TF1 *LandauGaussDist = new TF1("LandauGaussDist","[0]*exp((-0.5*((x-[1])/[2])^2))/([2]*sqrt(2*TMath::Pi()))*TMath::Landau(x,[3],[4])",xmin,xmax);
  cout<<"Par. shower: "<<0.5*h_NBINS*h_maxvalue<<", "<<h_MEAN<<", "<<h_RMS<<", "<<h_MEAN<<", "<<h_RMS<<endl;
  LandauGaussDist->SetParameter(0,0.5*h_NBINS*h_maxvalue);
  LandauGaussDist->SetParameter(1,h_MEAN);
  LandauGaussDist->SetParameter(2,h_RMS);
  LandauGaussDist->SetParameter(3,h_MEAN);
  LandauGaussDist->SetParameter(4,h_RMS);

  h->Fit(LandauGaussDist,"LRq");
  G_mu = LandauGaussDist->GetParameter(1);
  G_sig = LandauGaussDist->GetParameter(2);
  L_mpv = LandauGaussDist->GetParameter(3);
  L_fwhm = LandauGaussDist->GetParameter(4);

  MPV = LandauGaussDist->GetMaximumX();

  double maxvalue = LandauGaussDist->Eval(MPV);
  cout<<"Max. Value: "<<maxvalue<<", at x= "<<MPV<<endl;
  double x1=0.;
  double x2=MPV;
  for(int i=0;i<10000000;i++){
    x1 += 0.00001;
    double value = LandauGaussDist->Eval(x1);
    double distance = abs(maxvalue/2 - value);
    //cout<<distance<<endl;                                                                                                                                                                                                                                                                   
    if(distance < 0.01) break;
  }
  for(int i=0;i<10000000;i++){
    x2 += 0.00001;
    double value = LandauGaussDist->Eval(x2);
    double distance = abs(maxvalue/2 - value);
    //cout<<distance<<endl;                                                                                                                                                                                                                                                                   
    if(abs(maxvalue/2-value) < 0.01) break;
  }

  //cout<<"x1: "<<x1<<", x2: "<<x2<<endl;
  FWHM = x2 - x1;
  return;
}

void graph_setup_add(TGraph *g, string title, Color_t color){
    g->SetTitle(title.c_str());
    g->SetLineColor(color);
    g->SetLineWidth(3);
    return;
}

void analysis (string particle) {
    
    // double energies[N_ENERGIES] = {1, 2, 5, 8, 10, 20, 40, 60, 80, 100, 150};
    //double test_e[1]={20};
    //double test_e[2]={6.,80.};
    double test_e[N_ENERGIES] = {2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200};
    TVectorD energies(N_ENERGIES, test_e);
    TVectorD energies_tr(N_ENERGIES);
    for (int j = 0; j < N_ENERGIES; j++) energies_tr[j] = 1/TMath::Sqrt(energies[j]);
    // For conf11
    double W[15] = {4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6}; 
    TVectorD W_thicknesses(15, W);
    //Si: 650 650 650 650 500 500 500 500 500 500 320 320 320 320 320
    
    string filenames[N_ENERGIES];
    string base_path = "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/submit_jobs/LCIO2build_folder/LCIO2build_output/";
    // for (int j = 0; j < N_ENERGIES; j++) filenames[j] = base_path + "/CONF11/build/ECAL_QGSP_BERT_conf6_e-_" + to_string((int)round(energies[j])) +  "GeV_5kevt_build_masked.root";
    for (int j = 0; j < N_ENERGIES; j++) filenames[j] = base_path + "output_LCIO2Build_TB2022-06_"+particle+"_" + to_string((int)round(energies[j])) +  "GeV.root";
    
    // for (int j = 0; j < N_ENERGIES; j++) filenames[j] = base_path + "CONF11/build/ECAL_QGSP_BERT_conf8_e-_" + to_string((int)round(energies[j])) +  "GeV_5kevt_build_masked.root";
    //double test_zeros[1]={0.};
    //double test_zeros[2]={0.,0.};
    double test_zeros[N_ENERGIES] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // output_LCIO2Build_e_200GeV.root
    TVectorD zeros(N_ENERGIES, test_zeros);
    
    TVectorD mu_nhit(N_ENERGIES), mu_nhit_masked(N_ENERGIES), mu_sume(N_ENERGIES), mu_sume_masked(N_ENERGIES), mu_weight(N_ENERGIES), mu_weight_masked(N_ENERGIES);
    TVectorD sig_nhit(N_ENERGIES), sig_nhit_masked(N_ENERGIES), sig_sume(N_ENERGIES), sig_sume_masked(N_ENERGIES), sig_weight(N_ENERGIES), sig_weight_masked(N_ENERGIES);
    TVectorD res_nhit(N_ENERGIES), res_nhit_masked(N_ENERGIES), res_sume(N_ENERGIES), res_sume_masked(N_ENERGIES), res_weight(N_ENERGIES), res_weight_masked(N_ENERGIES);

    TVectorD mu_mol(N_ENERGIES), mu_mol_masked(N_ENERGIES);
    TVectorD sig_mol(N_ENERGIES), sig_mol_masked(N_ENERGIES);

    TVectorD mu_nhit_layer_0(N_ENERGIES), sig_nhit_layer_0(N_ENERGIES);
    TVectorD mu_nhit_layer_1(N_ENERGIES), sig_nhit_layer_1(N_ENERGIES);
    TVectorD mu_nhit_layer_2(N_ENERGIES), sig_nhit_layer_2(N_ENERGIES);
    TVectorD mu_nhit_layer_3(N_ENERGIES), sig_nhit_layer_3(N_ENERGIES);
    TVectorD mu_nhit_layer_4(N_ENERGIES), sig_nhit_layer_4(N_ENERGIES);
    TVectorD mu_nhit_layer_5(N_ENERGIES), sig_nhit_layer_5(N_ENERGIES);
    TVectorD mu_nhit_layer_6(N_ENERGIES), sig_nhit_layer_6(N_ENERGIES);
    TVectorD mu_nhit_layer_7(N_ENERGIES), sig_nhit_layer_7(N_ENERGIES);
    TVectorD mu_nhit_layer_8(N_ENERGIES), sig_nhit_layer_8(N_ENERGIES);
    TVectorD mu_nhit_layer_9(N_ENERGIES), sig_nhit_layer_9(N_ENERGIES);
    TVectorD mu_nhit_layer_10(N_ENERGIES), sig_nhit_layer_10(N_ENERGIES);
    TVectorD mu_nhit_layer_11(N_ENERGIES), sig_nhit_layer_11(N_ENERGIES);
    TVectorD mu_nhit_layer_12(N_ENERGIES), sig_nhit_layer_12(N_ENERGIES);
    TVectorD mu_nhit_layer_13(N_ENERGIES), sig_nhit_layer_13(N_ENERGIES);
    TVectorD mu_nhit_layer_14(N_ENERGIES), sig_nhit_layer_14(N_ENERGIES);
  
    TVectorD shower_nhit_max_layer(N_ENERGIES), shower_nhit_start_layer(N_ENERGIES), shower_nhit_end_layer(N_ENERGIES), shower_nhit_average(N_ENERGIES), shower_nhit_max(N_ENERGIES);

    TVectorD mu_nhit_layer_w_0(N_ENERGIES), sig_nhit_layer_w_0(N_ENERGIES);
    TVectorD mu_nhit_layer_w_1(N_ENERGIES), sig_nhit_layer_w_1(N_ENERGIES);
    TVectorD mu_nhit_layer_w_2(N_ENERGIES), sig_nhit_layer_w_2(N_ENERGIES);
    TVectorD mu_nhit_layer_w_3(N_ENERGIES), sig_nhit_layer_w_3(N_ENERGIES);
    TVectorD mu_nhit_layer_w_4(N_ENERGIES), sig_nhit_layer_w_4(N_ENERGIES);
    TVectorD mu_nhit_layer_w_5(N_ENERGIES), sig_nhit_layer_w_5(N_ENERGIES);
    TVectorD mu_nhit_layer_w_6(N_ENERGIES), sig_nhit_layer_w_6(N_ENERGIES);
    TVectorD mu_nhit_layer_w_7(N_ENERGIES), sig_nhit_layer_w_7(N_ENERGIES);
    TVectorD mu_nhit_layer_w_8(N_ENERGIES), sig_nhit_layer_w_8(N_ENERGIES);
    TVectorD mu_nhit_layer_w_9(N_ENERGIES), sig_nhit_layer_w_9(N_ENERGIES);
    TVectorD mu_nhit_layer_w_10(N_ENERGIES), sig_nhit_layer_w_10(N_ENERGIES);
    TVectorD mu_nhit_layer_w_11(N_ENERGIES), sig_nhit_layer_w_11(N_ENERGIES);
    TVectorD mu_nhit_layer_w_12(N_ENERGIES), sig_nhit_layer_w_12(N_ENERGIES);
    TVectorD mu_nhit_layer_w_13(N_ENERGIES), sig_nhit_layer_w_13(N_ENERGIES);
    TVectorD mu_nhit_layer_w_14(N_ENERGIES), sig_nhit_layer_w_14(N_ENERGIES);

  
    // Output filename
    TString result_name = "resolution_"+particle+"_result.root" ;
    TFile f(result_name, "recreate");

    for (int i_energy = 0; i_energy < N_ENERGIES; i_energy++) {
        cout << "Energy = " << energies[i_energy] << " GeV" << endl;
        string filename = filenames[i_energy];
        TFile *file = new TFile(filename.c_str(), "read");
        TTree *tree = (TTree*)file->Get("ecal");
    
        Long64_t nentries = tree->GetEntries();

	cout << "nentries: " << nentries <<endl;

        vector<float> *hit_energy = 0;
        vector<float> *hit_x = 0;
        vector<float> *hit_y = 0;
        vector<float> *hit_z = 0;
        vector<int> *hit_isMasked = 0;
        vector<int> *hit_slab = 0;
        
        TBranch *bhit_energy = 0;
        TBranch *bhit_x = 0;
        TBranch *bhit_y = 0;
        TBranch *bhit_z = 0;
        TBranch *bhit_isMasked = 0;
        TBranch *bhit_slab= 0;

        tree->SetBranchAddress("hit_energy", &hit_energy, &bhit_energy);
        tree->SetBranchAddress("hit_x", &hit_x, &bhit_x);
        tree->SetBranchAddress("hit_y", &hit_y, &bhit_y);
        tree->SetBranchAddress("hit_z", &hit_z, &bhit_z);
        tree->SetBranchAddress("hit_isMasked", &hit_isMasked, &bhit_isMasked);
        tree->SetBranchAddress("hit_slab", &hit_slab, &bhit_slab);

        // Resolution histos
        string e_str = to_string((int)round(energies[i_energy])) + "GeV";
	string part_string = particle;
	
        TH1I *h_nhit = new TH1I(("NumHits_" + part_string + "_" + e_str).c_str(), ("Number of hits " + part_string + " " + e_str).c_str(), 4000, 0, 4000);
        TH1I *h_nhit_masked = new TH1I(("NumHitsMask_" + part_string + "_" + e_str).c_str(), ("Number of hits, masked " + part_string + " " + e_str).c_str(), 4000, 0, 4000);
        TH1F *h_sume = new TH1F(("SumEnergy_" + part_string + "_" + e_str).c_str(), ("Sum Energy " + part_string + " " + e_str).c_str(), 500, 0, 15000);
        TH1F *h_sume_masked = new TH1F(("SumEnergyMask_" + part_string + "_" + e_str).c_str(), ("Sum Energy, masked " + part_string + " " + e_str).c_str(), 500, 0, 15000);
        TH1F *h_weight = new TH1F(("WSumEnergy_" + part_string + "_" + e_str).c_str(), ("W Sum Energy " + part_string + " " + e_str).c_str(), 500, 0, 70000);
        TH1F *h_weight_masked = new TH1F(("WSumEnergyMask_" + part_string + "_" + e_str).c_str(), ("W Sum Energy, masked " + part_string + " " + e_str).c_str(), 500, 0, 70000);
        
        // Moliere histos
        TH1F *h_mol = new TH1F(("Radius90_" + part_string + "_" + e_str).c_str(), ("Radius containing 90% of energy " + part_string + " " + e_str).c_str(), 100, 0., 100.);
        TH1F *h_mol_masked = new TH1F(("Radius90Masked_" + part_string + "_" + e_str).c_str(), ("Radius containing 90% of energy Masked " + part_string + " " + e_str).c_str(), 100, 0., 100.);

	// Shower profile histos
	TH1I *h_nhit_layer_0 = new TH1I(("NumHits_layer_0_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 0) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_1 = new TH1I(("NumHits_layer_1_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 1) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_2 = new TH1I(("NumHits_layer_2_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 2) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_3 = new TH1I(("NumHits_layer_3_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 3) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_4 = new TH1I(("NumHits_layer_4_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 4) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_5 = new TH1I(("NumHits_layer_5_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 5) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_6 = new TH1I(("NumHits_layer_6_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 6) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_7 = new TH1I(("NumHits_layer_7_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 7) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_8 = new TH1I(("NumHits_layer_8_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 8) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_9 = new TH1I(("NumHits_layer_9_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 9) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_10 = new TH1I(("NumHits_layer_10_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 10) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_11 = new TH1I(("NumHits_layer_11_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 11) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_12 = new TH1I(("NumHits_layer_12_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 12) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_13 = new TH1I(("NumHits_layer_13_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 13) " + part_string + " " + e_str).c_str(), 200, 0, 200);
	TH1I *h_nhit_layer_14 = new TH1I(("NumHits_layer_14_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 14) " + part_string + " " + e_str).c_str(), 200, 0, 200);

	TH1I *h_shower_nhit_max_layer = new TH1I(("ShowerNhitMaxLayer_" + part_string + "_" + e_str).c_str(), ("Shower Nhit Max. (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
	TH1I *h_shower_nhit_start_layer = new TH1I(("ShowerNhitStartLayer_" + part_string + "_" + e_str).c_str(), ("Shower Nhit Start (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
	TH1I *h_shower_nhit_end_layer = new TH1I(("ShowerNhitEndLayer_" + part_string + "_" + e_str).c_str(), ("Shower Nhit End (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
	TH1F *h_shower_nhit_average= new TH1F(("ShowerNhitAverage_" + part_string + "_" + e_str).c_str(), ("Shower Nhit Average  " + part_string + " " + e_str).c_str(), 200, 0, 200);
        TH1F *h_shower_nhit_max= new TH1F(("ShowerNhitMax_" + part_string + "_" + e_str).c_str(), ("Shower Nhit Max  " + part_string + " " + e_str).c_str(), 200, 0, 200);

	// weighted
	TH1I *h_nhit_layer_w_0 = new TH1I(("NumHits_layer_w_0_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 0) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_1 = new TH1I(("NumHits_layer_w_1_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 1) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_2 = new TH1I(("NumHits_layer_w_2_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 2) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_3 = new TH1I(("NumHits_layer_w_3_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 3) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_4 = new TH1I(("NumHits_layer_w_4_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 4) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_5 = new TH1I(("NumHits_layer_w_5_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 5) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_6 = new TH1I(("NumHits_layer_w_6_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 6) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_7 = new TH1I(("NumHits_layer_w_7_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 7) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_8 = new TH1I(("NumHits_layer_w_8_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 8) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_9 = new TH1I(("NumHits_layer_w_9_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 9) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_10 = new TH1I(("NumHits_layer_w_10_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 10) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_11 = new TH1I(("NumHits_layer_w_11_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 11) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_12 = new TH1I(("NumHits_layer_w_12_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 12) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_13 = new TH1I(("NumHits_layer_w_13_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 13) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1I *h_nhit_layer_w_14 = new TH1I(("NumHits_layer_w_14_" + part_string + "_" + e_str).c_str(), ("Number of hits weighted (layer 14) " + part_string + " " + e_str).c_str(), 100, 0, 1);

	//cout<<"debug 1"<<endl;
        for (int i_event = 0; i_event < (int)nentries; i_event++) {
	  // 10 for nentries
            // Resolution
            int nhit = 0;      int nhit_masked = 0;
            float sume = 0;    float sume_masked = 0;
            float weight = 0;  float weight_masked = 0;
	    //cout<<"debug 2"<<endl;
            tree->GetEntry(i_event);
	    if(i_event % 1000 == 0) cout << "Event " << to_string(i_event) << endl;
            //cout << "Event " << to_string(i_event) << endl;
	    //if (i_event == 3317) continue;
            
	    //cout<<"debug 3"<<endl;

	    get_res(nhit, nhit_masked, sume, sume_masked,
                    weight, weight_masked,
                    hit_energy, hit_slab, hit_isMasked,
                    W_thicknesses);

	    //cout<<"debug 4"<<endl;
	    
            h_nhit->Fill(nhit);        
	    h_nhit_masked->Fill(nhit_masked);
            h_sume->Fill(sume);         
	    h_sume_masked->Fill(sume_masked);
            h_weight->Fill(weight);     
	    h_weight_masked->Fill(weight_masked);

	    //cout<<"debug before fill mol"<<endl;
	    
            // Fill Moliere radii histograms
            h_mol->Fill(moliere(hit_energy, hit_x, hit_y, hit_z, hit_isMasked));
            h_mol_masked->Fill(moliere(hit_energy, hit_x, hit_y, hit_z, hit_isMasked, true));
            
	    //cout<<"debug after fill mol"<<endl;
	    
	    // Fill shower profile
	    h_nhit_layer_0->Fill(hit_layer(0, hit_slab, hit_isMasked));
	    h_nhit_layer_1->Fill(hit_layer(1, hit_slab, hit_isMasked));
	    h_nhit_layer_2->Fill(hit_layer(2, hit_slab, hit_isMasked));
	    h_nhit_layer_3->Fill(hit_layer(3, hit_slab, hit_isMasked));
	    h_nhit_layer_4->Fill(hit_layer(4, hit_slab, hit_isMasked));
	    h_nhit_layer_5->Fill(hit_layer(5, hit_slab, hit_isMasked));
	    h_nhit_layer_6->Fill(hit_layer(6, hit_slab, hit_isMasked));
	    h_nhit_layer_7->Fill(hit_layer(7, hit_slab, hit_isMasked));
	    h_nhit_layer_8->Fill(hit_layer(8, hit_slab, hit_isMasked));
	    h_nhit_layer_9->Fill(hit_layer(9, hit_slab, hit_isMasked));
	    h_nhit_layer_10->Fill(hit_layer(10, hit_slab, hit_isMasked));
	    h_nhit_layer_11->Fill(hit_layer(11, hit_slab, hit_isMasked));
	    h_nhit_layer_12->Fill(hit_layer(12, hit_slab, hit_isMasked));
	    h_nhit_layer_13->Fill(hit_layer(13, hit_slab, hit_isMasked));
	    h_nhit_layer_14->Fill(hit_layer(14, hit_slab, hit_isMasked));
	    
	    h_nhit_layer_w_0->Fill(hit_layer(0, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_1->Fill(hit_layer(1, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_2->Fill(hit_layer(2, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_3->Fill(hit_layer(3, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_4->Fill(hit_layer(4, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_5->Fill(hit_layer(5, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_6->Fill(hit_layer(6, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_7->Fill(hit_layer(7, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_8->Fill(hit_layer(8, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_9->Fill(hit_layer(9, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_10->Fill(hit_layer(10, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_11->Fill(hit_layer(11, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_12->Fill(hit_layer(12, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_13->Fill(hit_layer(13, hit_slab, hit_isMasked, false, true));
            h_nhit_layer_w_14->Fill(hit_layer(14, hit_slab, hit_isMasked, false, true));

	    double nhit_shower_maxvalue = 0.;
	    int nhit_ilayermax = -1;
	    int nhit_ilayerstart = -1;
	    int nhit_ilayerend = -1;
	    double nhit_shower_maxvalue_w = 0.;

	    for(int ilayer=0; ilayer<15; ilayer++){
              double thislayer = hit_layer(ilayer, hit_slab, hit_isMasked);
              if(thislayer > 0.){
                nhit_ilayerstart = ilayer;
                break;
              }
            }
	    
	    for(int ilayer=0; ilayer<15; ilayer++){
	      double thislayer = hit_layer(ilayer, hit_slab, hit_isMasked);
	      double thislayer_w = hit_layer(ilayer, hit_slab, hit_isMasked, false, true);
	      if(thislayer > 0.) {
		if(thislayer > nhit_shower_maxvalue){
		  nhit_shower_maxvalue = thislayer;
		  nhit_shower_maxvalue_w = thislayer_w;
		  nhit_ilayermax = ilayer;
		}
	      }
	    }
	    
	    double layertempmax = nhit_shower_maxvalue;
	    for(int ilayer=nhit_ilayermax; ilayer<15; ilayer++){
              if(nhit_shower_maxvalue > 5.) {
		double thislayer = hit_layer(ilayer, hit_slab, hit_isMasked);
		double drop = (layertempmax - thislayer)/layertempmax; 
		if( (drop > 0.) and (layertempmax > 1.)) layertempmax = thislayer;
		else {
		  nhit_ilayerend = ilayer - 1; 
		}  
	      }
	    } 

	    double nhit_showersum = 0;
	    for(int ilayer=0; ilayer<15; ilayer++){
              double thislayer = hit_layer(ilayer, hit_slab, hit_isMasked);
	      nhit_showersum += thislayer;
	    }
	    double nhit_shower_averagevalue = nhit_showersum/15;
	    	    	    
	    h_shower_nhit_max_layer->Fill(nhit_ilayermax);
	    h_shower_nhit_start_layer->Fill(nhit_ilayerstart);
	    h_shower_nhit_end_layer->Fill(nhit_ilayerend);
	    h_shower_nhit_average->Fill(nhit_shower_averagevalue);
	    h_shower_nhit_max->Fill(nhit_shower_maxvalue);
	    
            hit_isMasked->clear();
            hit_energy->clear();
            hit_slab->clear();
	
        }

        // Resolution
        fit_res(h_nhit, mu_nhit[i_energy], sig_nhit[i_energy], res_nhit[i_energy]);
        fit_res(h_nhit_masked, mu_nhit_masked[i_energy], sig_nhit_masked[i_energy], res_nhit_masked[i_energy]);
        fit_res(h_sume, mu_sume[i_energy], sig_sume[i_energy], res_sume[i_energy]);
        fit_res(h_sume_masked, mu_sume_masked[i_energy], sig_sume_masked[i_energy], res_sume_masked[i_energy]);
        fit_res(h_weight, mu_weight[i_energy], sig_weight[i_energy], res_weight[i_energy]);
        fit_res(h_weight_masked, mu_weight_masked[i_energy], sig_weight_masked[i_energy], res_weight_masked[i_energy]);

	// Moliere fits
	//NEW method
	// g_mu, g_sig, l_mpv, l_fwhm
	double parmolfit[4] = {20.,40.,50.,10.};
	double parmolfit_masked[4] = {20.,40.,50.,10.};
	// MPV and FWHM
	double result_molfit[2] = {0.,0.};
	double result_molfit_masked[2] = {0.,0.};
	// moliere min and moliere max range
	double range[2] = {0.1, 100.};
	
	//cout<<"before mol fit"<<endl;
	//fit_mol(h_mol, mu, sig, mpv, fwhm, range[0], range[1]);
	//h_mol->Fit("gaus");
	fit_mol(h_mol, parmolfit[0], parmolfit[1], parmolfit[2], parmolfit[3], result_molfit[0], result_molfit[1], range[0], range[1]);
	mu_mol[i_energy] = result_molfit[0];
	sig_mol[i_energy] = 0.5*result_molfit[1];
	
	fit_mol(h_mol_masked, parmolfit_masked[0], parmolfit_masked[1], parmolfit_masked[2], parmolfit_masked[3], result_molfit_masked[0], result_molfit_masked[1], range[0], range[1]);
        mu_mol_masked[i_energy] = result_molfit_masked[0];
        sig_mol_masked[i_energy] = 0.5*result_molfit_masked[1];

	cout<<"Molière radius: ("<<mu_mol[i_energy]<<" +- "<<sig_mol[i_energy]<<") mm"<<endl;
	cout<<"Molière radius (masked): ("<<mu_mol_masked[i_energy]<<" +- "<<sig_mol_masked[i_energy]<<") mm"<<endl;
	
	// Shower profile fits
	// Same method than moliere
	// g_mu, g_sig, l_mpv, l_fwhm
        double parshofit[4] = {4.,2.,2.,2.};
        double result_shofit[15][2];
	double rangesho[2] = {0., 200.};
        fit_shower(h_nhit_layer_0, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[0][0], result_shofit[0][1], rangesho[0], rangesho[1]);
        fit_shower(h_nhit_layer_1, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[1][0], result_shofit[1][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_2, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[2][0], result_shofit[2][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_3, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[3][0], result_shofit[3][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_4, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[4][0], result_shofit[4][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_5, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[5][0], result_shofit[5][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_6, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[6][0], result_shofit[6][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_7, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[7][0], result_shofit[7][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_8, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[8][0], result_shofit[8][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_9, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[9][0], result_shofit[9][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_10, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[10][0], result_shofit[10][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_11, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[11][0], result_shofit[11][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_12, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[12][0], result_shofit[12][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_13, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[13][0], result_shofit[13][1], rangesho[0], rangesho[1]);
	fit_shower(h_nhit_layer_14, parshofit[0], parshofit[1], parshofit[2], parshofit[3], result_shofit[14][0], result_shofit[14][1], rangesho[0], rangesho[1]);
	
	mu_nhit_layer_0[i_energy] = result_shofit[0][0];
        sig_nhit_layer_0[i_energy] = 0.5*result_shofit[0][1];
	mu_nhit_layer_1[i_energy] = result_shofit[1][0];
        sig_nhit_layer_1[i_energy] = 0.5*result_shofit[1][1];
	mu_nhit_layer_2[i_energy] = result_shofit[2][0];
        sig_nhit_layer_2[i_energy] = 0.5*result_shofit[2][1];
	mu_nhit_layer_3[i_energy] = result_shofit[3][0];
        sig_nhit_layer_3[i_energy] = 0.5*result_shofit[3][1];
	mu_nhit_layer_4[i_energy] = result_shofit[4][0];
        sig_nhit_layer_4[i_energy] = 0.5*result_shofit[4][1];
	mu_nhit_layer_5[i_energy] = result_shofit[5][0];
        sig_nhit_layer_5[i_energy] = 0.5*result_shofit[5][1];
	mu_nhit_layer_6[i_energy] = result_shofit[6][0];
        sig_nhit_layer_6[i_energy] = 0.5*result_shofit[6][1];
	mu_nhit_layer_7[i_energy] = result_shofit[7][0];
        sig_nhit_layer_7[i_energy] = 0.5*result_shofit[7][1];
	mu_nhit_layer_8[i_energy] = result_shofit[8][0];
        sig_nhit_layer_8[i_energy] = 0.5*result_shofit[8][1];
	mu_nhit_layer_9[i_energy] = result_shofit[9][0];
        sig_nhit_layer_9[i_energy] = 0.5*result_shofit[9][1];
	mu_nhit_layer_10[i_energy] = result_shofit[10][0];
        sig_nhit_layer_10[i_energy] = 0.5*result_shofit[10][1];
	mu_nhit_layer_11[i_energy] = result_shofit[11][0];
        sig_nhit_layer_11[i_energy] = 0.5*result_shofit[11][1];
	mu_nhit_layer_12[i_energy] = result_shofit[12][0];
        sig_nhit_layer_12[i_energy] = 0.5*result_shofit[12][1];
	mu_nhit_layer_13[i_energy] = result_shofit[13][0];
        sig_nhit_layer_13[i_energy] = 0.5*result_shofit[13][1];
	mu_nhit_layer_14[i_energy] = result_shofit[14][0];
        sig_nhit_layer_14[i_energy] = 0.5*result_shofit[14][1];
    
	shower_nhit_max_layer[i_energy] = h_shower_nhit_max_layer->GetMean();
	shower_nhit_start_layer[i_energy] = h_shower_nhit_start_layer->GetMean();
	shower_nhit_end_layer[i_energy] = h_shower_nhit_end_layer->GetMean(); //-2 to include layer -1 and because hist starts in 1
	shower_nhit_average[i_energy] = h_shower_nhit_average->GetMean(); //-1 because hists starts in 1
	shower_nhit_max[i_energy] = h_shower_nhit_max->GetMean();
	
	cout<<"Shower layer start: "<<shower_nhit_start_layer[i_energy]<<", max: "<<shower_nhit_max_layer[i_energy]<<", end: "<<shower_nhit_end_layer[i_energy]<<endl;
	cout<<"Shower layer average: "<<shower_nhit_average[i_energy]<<", max. value: "<<shower_nhit_max[i_energy]<<endl;;
	  
	// Shower profile fits
        // Same method than moliere
        // g_mu, g_sig, l_mpv, l_fwhm                                                                                                                                                                                     
        double parshofit_w[4] = {1.,1.,1.,1.};
        double result_shofit_w[15][2];
        double rangesho_w[2] = {0., 200.};
        fit_shower(h_nhit_layer_w_0, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[0][0], result_shofit_w[0][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_1, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[1][0], result_shofit_w[1][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_2, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[2][0], result_shofit_w[2][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_3, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[3][0], result_shofit_w[3][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_4, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[4][0], result_shofit_w[4][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_5, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[5][0], result_shofit_w[5][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_6, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[6][0], result_shofit_w[6][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_7, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[7][0], result_shofit_w[7][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_8, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[8][0], result_shofit_w[8][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_9, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[9][0], result_shofit_w[9][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_10, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[10][0], result_shofit_w[10][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_11, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[11][0], result_shofit_w[11][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_12, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[12][0], result_shofit_w[12][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_13, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[13][0], result_shofit_w[13][1], rangesho_w[0], rangesho_w[1]);
        fit_shower(h_nhit_layer_w_14, parshofit_w[0], parshofit_w[1], parshofit_w[2], parshofit_w[3], result_shofit_w[14][0], result_shofit_w[14][1], rangesho_w[0], rangesho_w[1]);

        mu_nhit_layer_w_0[i_energy] = result_shofit_w[0][0];
        sig_nhit_layer_w_0[i_energy] = 0.5*result_shofit_w[0][1];
        mu_nhit_layer_w_1[i_energy] = result_shofit_w[1][0];
        sig_nhit_layer_w_1[i_energy] = 0.5*result_shofit_w[1][1];
        mu_nhit_layer_w_2[i_energy] = result_shofit_w[2][0];
        sig_nhit_layer_w_2[i_energy] = 0.5*result_shofit_w[2][1];
        mu_nhit_layer_w_3[i_energy] = result_shofit_w[3][0];
        sig_nhit_layer_w_3[i_energy] = 0.5*result_shofit_w[3][1];
        mu_nhit_layer_w_4[i_energy] = result_shofit_w[4][0];
        sig_nhit_layer_w_4[i_energy] = 0.5*result_shofit_w[4][1];
        mu_nhit_layer_w_5[i_energy] = result_shofit_w[5][0];
        sig_nhit_layer_w_5[i_energy] = 0.5*result_shofit_w[5][1];
        mu_nhit_layer_w_6[i_energy] = result_shofit_w[6][0];
        sig_nhit_layer_w_6[i_energy] = 0.5*result_shofit_w[6][1];
        mu_nhit_layer_w_7[i_energy] = result_shofit_w[7][0];
        sig_nhit_layer_w_7[i_energy] = 0.5*result_shofit_w[7][1];
        mu_nhit_layer_w_8[i_energy] = result_shofit_w[8][0];
        sig_nhit_layer_w_8[i_energy] = 0.5*result_shofit_w[8][1];
        mu_nhit_layer_w_9[i_energy] = result_shofit_w[9][0];
        sig_nhit_layer_w_9[i_energy] = 0.5*result_shofit_w[9][1];
        mu_nhit_layer_w_10[i_energy] = result_shofit_w[10][0];
        sig_nhit_layer_w_10[i_energy] = 0.5*result_shofit_w[10][1];
        mu_nhit_layer_w_11[i_energy] = result_shofit_w[11][0];
        sig_nhit_layer_w_11[i_energy] = 0.5*result_shofit_w[11][1];
        mu_nhit_layer_w_12[i_energy] = result_shofit_w[12][0];
        sig_nhit_layer_w_12[i_energy] = 0.5*result_shofit_w[12][1];
        mu_nhit_layer_w_13[i_energy] = result_shofit_w[13][0];
        sig_nhit_layer_w_13[i_energy] = 0.5*result_shofit_w[13][1];
        mu_nhit_layer_w_14[i_energy] = result_shofit_w[14][0];
        sig_nhit_layer_w_14[i_energy] = 0.5*result_shofit_w[14][1];


	//Write objects
        f.WriteTObject(h_nhit);    f.WriteTObject(h_nhit_masked);
        f.WriteTObject(h_sume);    f.WriteTObject(h_sume_masked);
        f.WriteTObject(h_weight);  f.WriteTObject(h_weight_masked);
        f.WriteTObject(h_mol);     f.WriteTObject(h_mol_masked);
	
	f.WriteTObject(h_nhit_layer_0);
	f.WriteTObject(h_nhit_layer_1);
	f.WriteTObject(h_nhit_layer_2);
	f.WriteTObject(h_nhit_layer_3);
	f.WriteTObject(h_nhit_layer_4);
	f.WriteTObject(h_nhit_layer_5);
	f.WriteTObject(h_nhit_layer_6);
	f.WriteTObject(h_nhit_layer_7);
	f.WriteTObject(h_nhit_layer_8);
	f.WriteTObject(h_nhit_layer_9);
	f.WriteTObject(h_nhit_layer_10);
	f.WriteTObject(h_nhit_layer_11);
	f.WriteTObject(h_nhit_layer_12);
	f.WriteTObject(h_nhit_layer_13);
	f.WriteTObject(h_nhit_layer_14);
	f.WriteTObject(h_shower_nhit_max_layer);
	f.WriteTObject(h_shower_nhit_start_layer);
	f.WriteTObject(h_shower_nhit_end_layer);
	f.WriteTObject(h_shower_nhit_average);
	f.WriteTObject(h_shower_nhit_max);
	
	f.WriteTObject(h_nhit_layer_w_0);
        f.WriteTObject(h_nhit_layer_w_1);
        f.WriteTObject(h_nhit_layer_w_2);
        f.WriteTObject(h_nhit_layer_w_3);
        f.WriteTObject(h_nhit_layer_w_4);
        f.WriteTObject(h_nhit_layer_w_5);
        f.WriteTObject(h_nhit_layer_w_6);
        f.WriteTObject(h_nhit_layer_w_7);
        f.WriteTObject(h_nhit_layer_w_8);
        f.WriteTObject(h_nhit_layer_w_9);
        f.WriteTObject(h_nhit_layer_w_10);
        f.WriteTObject(h_nhit_layer_w_11);
        f.WriteTObject(h_nhit_layer_w_12);
        f.WriteTObject(h_nhit_layer_w_13);
        f.WriteTObject(h_nhit_layer_w_14);
    }
    
    f.WriteObject(&mu_nhit, "mu_nhit");                     f.WriteObject(&sig_nhit, "sig_nhit");                   f.WriteObject(&res_nhit, "res_nhit");
    f.WriteObject(&mu_nhit_masked, "mu_nhit_masked");       f.WriteObject(&sig_nhit_masked, "sig_nhit_masked");     f.WriteObject(&res_nhit_masked, "res_nhit_masked");
    f.WriteObject(&mu_sume, "mu_sume");                     f.WriteObject(&sig_sume, "sig_sume");                   f.WriteObject(&res_sume, "res_sume");
    f.WriteObject(&mu_sume_masked, "mu_sume_masked");       f.WriteObject(&sig_sume_masked, "sig_sume_masked");     f.WriteObject(&res_sume_masked, "res_sume_masked");
    f.WriteObject(&mu_weight, "mu_weight");                 f.WriteObject(&sig_weight, "sig_weight");               f.WriteObject(&res_weight, "res_weight");
    f.WriteObject(&mu_weight_masked, "mu_weight_masked");   f.WriteObject(&sig_weight_masked, "sig_weight_masked"); f.WriteObject(&res_weight_masked, "res_weight_masked");
    f.WriteObject(&mu_mol, "mu_mol");   f.WriteObject(&mu_mol_masked, "mu_mol_masked"); f.WriteObject(&sig_mol, "sig_mol");   f.WriteObject(&sig_mol_masked, "sig_mol_masked");
    
    f.WriteTObject(&mu_nhit_layer_0, "mu_nhit_layer_0");
    f.WriteTObject(&mu_nhit_layer_1, "mu_nhit_layer_1");
    f.WriteTObject(&mu_nhit_layer_2, "mu_nhit_layer_2");
    f.WriteTObject(&mu_nhit_layer_3, "mu_nhit_layer_3");
    f.WriteTObject(&mu_nhit_layer_4, "mu_nhit_layer_4");
    f.WriteTObject(&mu_nhit_layer_5, "mu_nhit_layer_5");
    f.WriteTObject(&mu_nhit_layer_6, "mu_nhit_layer_6");
    f.WriteTObject(&mu_nhit_layer_7, "mu_nhit_layer_7");
    f.WriteTObject(&mu_nhit_layer_8, "mu_nhit_layer_8");
    f.WriteTObject(&mu_nhit_layer_9, "mu_nhit_layer_9");
    f.WriteTObject(&mu_nhit_layer_10, "mu_nhit_layer_10");
    f.WriteTObject(&mu_nhit_layer_11, "mu_nhit_layer_11");
    f.WriteTObject(&mu_nhit_layer_12, "mu_nhit_layer_12");
    f.WriteTObject(&mu_nhit_layer_13, "mu_nhit_layer_13");
    f.WriteTObject(&mu_nhit_layer_14, "mu_nhit_layer_14");
    f.WriteTObject(&sig_nhit_layer_0, "sig_nhit_layer_0");
    f.WriteTObject(&sig_nhit_layer_1, "sig_nhit_layer_1");
    f.WriteTObject(&sig_nhit_layer_2, "sig_nhit_layer_2");
    f.WriteTObject(&sig_nhit_layer_3, "sig_nhit_layer_3");
    f.WriteTObject(&sig_nhit_layer_4, "sig_nhit_layer_4");
    f.WriteTObject(&sig_nhit_layer_5, "sig_nhit_layer_5");
    f.WriteTObject(&sig_nhit_layer_6, "sig_nhit_layer_6");
    f.WriteTObject(&sig_nhit_layer_7, "sig_nhit_layer_7");
    f.WriteTObject(&sig_nhit_layer_8, "sig_nhit_layer_8");
    f.WriteTObject(&sig_nhit_layer_9, "sig_nhit_layer_9");
    f.WriteTObject(&sig_nhit_layer_10, "sig_nhit_layer_10");
    f.WriteTObject(&sig_nhit_layer_11, "sig_nhit_layer_11");
    f.WriteTObject(&sig_nhit_layer_12, "sig_nhit_layer_12");
    f.WriteTObject(&sig_nhit_layer_13, "sig_nhit_layer_13");
    f.WriteTObject(&sig_nhit_layer_14, "sig_nhit_layer_14");
    f.WriteTObject(&shower_nhit_max_layer, "shower_nhit_max_layer");
    f.WriteTObject(&shower_nhit_start_layer, "shower_nhit_start_layer");
    f.WriteTObject(&shower_nhit_end_layer, "shower_nhit_end_layer");
    f.WriteTObject(&shower_nhit_average, "shower_nhit_average");
    f.WriteTObject(&shower_nhit_max, "shower_nhit_max");

    f.WriteTObject(&mu_nhit_layer_w_0, "mu_nhit_layer_w_0");
    f.WriteTObject(&mu_nhit_layer_w_1, "mu_nhit_layer_w_1");
    f.WriteTObject(&mu_nhit_layer_w_2, "mu_nhit_layer_w_2");
    f.WriteTObject(&mu_nhit_layer_w_3, "mu_nhit_layer_w_3");
    f.WriteTObject(&mu_nhit_layer_w_4, "mu_nhit_layer_w_4");
    f.WriteTObject(&mu_nhit_layer_w_5, "mu_nhit_layer_w_5");
    f.WriteTObject(&mu_nhit_layer_w_6, "mu_nhit_layer_w_6");
    f.WriteTObject(&mu_nhit_layer_w_7, "mu_nhit_layer_w_7");
    f.WriteTObject(&mu_nhit_layer_w_8, "mu_nhit_layer_w_8");
    f.WriteTObject(&mu_nhit_layer_w_9, "mu_nhit_layer_w_9");
    f.WriteTObject(&mu_nhit_layer_w_10, "mu_nhit_layer_w_10");
    f.WriteTObject(&mu_nhit_layer_w_11, "mu_nhit_layer_w_11");
    f.WriteTObject(&mu_nhit_layer_w_12, "mu_nhit_layer_w_12");
    f.WriteTObject(&mu_nhit_layer_w_13, "mu_nhit_layer_w_13");
    f.WriteTObject(&mu_nhit_layer_w_14, "mu_nhit_layer_w_14");
    f.WriteTObject(&sig_nhit_layer_w_0, "sig_nhit_layer_w_0");
    f.WriteTObject(&sig_nhit_layer_w_1, "sig_nhit_layer_w_1");
    f.WriteTObject(&sig_nhit_layer_w_2, "sig_nhit_layer_w_2");
    f.WriteTObject(&sig_nhit_layer_w_3, "sig_nhit_layer_w_3");
    f.WriteTObject(&sig_nhit_layer_w_4, "sig_nhit_layer_w_4");
    f.WriteTObject(&sig_nhit_layer_w_5, "sig_nhit_layer_w_5");
    f.WriteTObject(&sig_nhit_layer_w_6, "sig_nhit_layer_w_6");
    f.WriteTObject(&sig_nhit_layer_w_7, "sig_nhit_layer_w_7");
    f.WriteTObject(&sig_nhit_layer_w_8, "sig_nhit_layer_w_8");
    f.WriteTObject(&sig_nhit_layer_w_9, "sig_nhit_layer_w_9");
    f.WriteTObject(&sig_nhit_layer_w_10, "sig_nhit_layer_w_10");
    f.WriteTObject(&sig_nhit_layer_w_11, "sig_nhit_layer_w_11");
    f.WriteTObject(&sig_nhit_layer_w_12, "sig_nhit_layer_w_12");
    f.WriteTObject(&sig_nhit_layer_w_13, "sig_nhit_layer_w_13");
    f.WriteTObject(&sig_nhit_layer_w_14, "sig_nhit_layer_w_14");
    
    f.WriteObject(&energies, "energies");
    f.WriteObject(&energies_tr, "energies_tr");
    f.WriteObject(&W_thicknesses, "W_thicknesses");
    f.Close();

}
