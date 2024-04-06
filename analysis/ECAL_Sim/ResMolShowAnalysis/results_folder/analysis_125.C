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

#define N_ENERGIES 1 //18
#define N_ECAL_LAYERS 15

void get_res(int &nhit, float &sume, float &weight, vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses, vector<int> *hit_isMasked, bool &masked) {
  if(hit_energy->size() > 0){    
    //cout<<W_thicknesses.Min()<<endl;
    // First option: use the minimum of the used
    // Second option: use the paper as reference 0.4X0, X0=3.5mm
    // It was weight_masked += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
    // New version ?
    
   for (int j = 0; j < hit_energy->size(); j++) {
     if( masked && hit_isMasked->at(j) == 1 ) continue;
            nhit += 1;
            sume += hit_energy->at(j);
            weight += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
   }
  }
    return;
}

void hits_layer(float hlv[N_ECAL_LAYERS], vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses, vector<int> *hit_isMasked, bool masked=false, bool normalized=false, string count_type="nhit"){
  float hit_count[N_ECAL_LAYERS];
  float sume[N_ECAL_LAYERS];
  float sume_w[N_ECAL_LAYERS];
  
  float sume_total = 1.;
  float sume_w_total = 1.;

  float weight = 1.;

  // Initial values
  for (int ilayer=0; ilayer < N_ECAL_LAYERS; ilayer++){
    hit_count[ilayer] = 0. ;
    sume[ilayer] = 0. ;
    sume_w[ilayer] = 0. ;
  }
  
  if(hit_energy->size() > 0){
    for (int j = 0; j < hit_energy->size(); j++) {
      if( masked && hit_isMasked->at(j) == 1 ) continue;
      sume_total += hit_energy->at(j);
      sume_w_total += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
      //cout<<"|DEBUG(hit_layer)| sume total: "<<sume_total<<", sume_w_total: "<<sume_w_total<<endl;
    }
  }
  
  for (int ilayer=0; ilayer < N_ECAL_LAYERS; ilayer++) {
    if(hit_slab->size() > 0){
      for( int j = 0; j < hit_energy->size(); j++ ) {
	if( masked && hit_isMasked->at(j) == 1 ) continue;
	if( hit_slab->at(j) == ilayer ) {
	  hit_count[ilayer] += 1;
	  sume[ilayer] += hit_energy->at(j);
	  sume_w[ilayer] += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
	}
      }
      if(normalized == true) {
	if(count_type == "nhit") weight = 1./hit_energy->size();
	if(count_type == "sume") weight = 1./sume_total;
	if(count_type == "weight") weight = 1./sume_w_total;
      }  
      
      if(count_type == "nhit") hlv[ilayer] = hit_count[ilayer]*weight;
      if(count_type == "sume") hlv[ilayer] = sume[ilayer]*weight;
      if(count_type == "weight") hlv[ilayer] = sume_w[ilayer]*weight;
    }
    else hlv[ilayer] = 0.;
  
  }
  //cout<<"| DEBUG(hits_layer) | "<<"Mode: "<<count_type<<", normalized: "<<normalized<<endl;
  //cout<<"| DEBUG(hits_layer) | "<<"Hit energy size: "<<hit_energy->size()<<". Hit slab size: "<<hit_slab->size()<<endl;
  //cout<<"| DEBUG(hits_layer) | "<<"hlv: "<<hlv[0]<<" "<<hlv[1]<<" "<<hlv[2]<<" "<<hlv[3]<<" "<<hlv[4]<<" "<<hlv[5]<<" "<<hlv[6]<<" "<<hlv[7]<<" "<<hlv[8]<<" "<<hlv[9]<<" "<<hlv[10]<<" "<<hlv[11]<<" "<<hlv[12]<<" "<<hlv[13]<<" "<<hlv[14]<<endl;
    
  
  return 0;
}

void shower_variables(float entries, float array[N_ECAL_LAYERS], float array_n[N_ECAL_LAYERS], float &shower_maxvalue, float &shower_maxvalue_n, int &ilayermax,
		      int &ilayerstart, int &ilayerstart_10, int &ilayerend, int &ilayerend_10, string count_type = "nhit") {
  
  float threshold = 0.;
  float percentage = 0.1;

  if(count_type == "nhit") threshold = 2.;
  if(count_type == "weight") threshold = 2.;  
  /*
  cout<<count_type<<" shower variables"<<endl;
  for(int ilayer=0; ilayer<N_ECAL_LAYERS; ilayer++)cout<<array[ilayer]<<" ";
  cout<<""<<endl;
  */
  if(entries > 0){
    for(int ilayer=0; ilayer<N_ECAL_LAYERS; ilayer++){
      float thislayer = array[ilayer];
      float thislayer_n = array_n[ilayer];
      if(thislayer > shower_maxvalue){
	shower_maxvalue = thislayer;
	shower_maxvalue_n = thislayer_n;
	ilayermax = ilayer;
      }
    }

    if((shower_maxvalue_n > 1/N_ECAL_LAYERS) && (shower_maxvalue > threshold)){
      for(int ilayer=N_ECAL_LAYERS-1; ilayer>ilayermax; ilayer--){
	float thislayer = array[ilayer];
	if((thislayer > threshold) && (thislayer < 0.5*shower_maxvalue)) {
	  ilayerend = ilayer;
	  break;
	}
      }
      for(int ilayer=N_ECAL_LAYERS-1; ilayer>ilayermax; ilayer--){
        float thislayer = array[ilayer];
	if((thislayer > threshold) && (thislayer > percentage*shower_maxvalue) && (thislayer < 0.3*shower_maxvalue)) {
	  ilayerend_10 = ilayer;
	  break;
	}
      }
      //cout<<count_type<<" end 0.1: "<<ilayerend_10<<", end: "<<ilayerend<<endl;

      for(int ilayer=0; ilayer<ilayermax; ilayer++){
	float thislayer = array[ilayer];
	//cout<<"|start loop| layer "<<ilayer<<" value: "<<thislayer<<endl;
	if((thislayer > threshold) && (thislayer < 0.5*shower_maxvalue)) {
	  ilayerstart = ilayer;
	  break;
	}
      }
      for(int ilayer=0; ilayer<ilayermax; ilayer++){
        float thislayer = array[ilayer];
	//cout<<"|start loop per| layer "<<ilayer<<" value: "<<thislayer<<endl;
        if((thislayer > threshold) && (thislayer > percentage*shower_maxvalue) && (thislayer < 0.3*shower_maxvalue)) {
	  ilayerstart_10 = ilayer;
	  break;
	}
      }
      //cout<<count_type<<" start 0.1: "<<ilayerstart_10<<", start: "<<ilayerstart<<endl;
    }
  }
  return 0;
}

float moliere(vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses,
              vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z,
              vector<int> * hit_isMasked, bool masked=false, float containment = 0.90) {

    vector<float> hit_rs;
    vector<float> hit_es;
    float mol_rad = 0.;
    float sume = 0.;
    float wx = 0.; float wy = 0.; float wz = 0.;
    float r = 0.;
    if(hit_energy->size() > 0){
      //cout<<"before first loop in molfunc"<<endl;
    
      for (int j = 0; j < hit_energy->size(); j++){
        if (masked && hit_isMasked->at(j) == 1) continue;
        sume += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
        wx += hit_x->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
        wy += hit_y->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
        wz += hit_z->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
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
      for (auto k: sort_indexes(hit_rs)) hit_es.push_back(hit_energy->at(k) * W_thicknesses[hit_slab->at(k)]/(0.4*3.5));
      
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

void barycenter(vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses, 
		vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z, float bar_xyzr[4], 
		vector<int> * hit_isMasked, bool masked=false) {
  
  float sume = 0.;
  float wx = 0.; float wy = 0.; float wz = 0.;
  float bary_x = 0., bary_y = 0., bary_z = 0.;
  if(hit_energy->size() > 0){
    for (int j = 0; j < hit_energy->size(); j++){
      if (masked && hit_isMasked->at(j) == 1) continue;
      //cout<<"DBG "<<j<<endl;
      sume += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
      wx += hit_x->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
      wy += hit_y->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
      wz += hit_z->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
    }
    
    bary_x = wx / sume;
    bary_y = wy / sume;
    bary_z = wz / sume;
  
    if( (sume == 0) or (sume < 0) ){
      bary_x = 0.;
      bary_y = 0.;
      bary_z = 0.;
    }
  }
  
  bar_xyzr[0] = bary_x;
  bar_xyzr[1] = bary_y;
  bar_xyzr[2] = bary_z;
  bar_xyzr[3] = pow(pow(bary_x,2)+pow(bary_y,2),0.5);

  return 0;
}

void bary_layer(float blv[N_ECAL_LAYERS][3], vector<float> * hit_energy, vector<int> *hit_slab, 
		TVectorD W_thicknesses, vector<float> * hit_x, vector<float> * hit_y, 
		vector<float> * hit_z, vector<int> * hit_isMasked, bool masked=false) {

  float sume_w[N_ECAL_LAYERS];
  float wx[N_ECAL_LAYERS]; float wy[N_ECAL_LAYERS];
  float bary_x[N_ECAL_LAYERS]; float bary_y[N_ECAL_LAYERS];
  
  for (int ilayer=0; ilayer < N_ECAL_LAYERS; ilayer++){
    sume_w[ilayer] = 0. ;
    wx[ilayer] = 0.;
    wy[ilayer] = 0.;
    bary_x[ilayer] = 0.;
    bary_y[ilayer] = 0.;
    blv[ilayer][0] = 0.;
    blv[ilayer][1] = 0.;
    blv[ilayer][2] = 0.;
  }
  
  if(hit_energy->size() > 0){
    for (int ilayer=0; ilayer < N_ECAL_LAYERS; ilayer++) {
      if(hit_slab->size() > 0){
	for (int j = 0; j < hit_energy->size(); j++) {
	  if( masked && hit_isMasked->at(j) == 1 ) continue;
	  if( hit_slab->at(j) == ilayer ) {
	    sume_w[ilayer] += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
	  }
	}
	for (int j = 0; j < hit_energy->size(); j++) {
	  if( masked && hit_isMasked->at(j) == 1 ) continue;
	  if( hit_slab->at(j) == ilayer ) { 
	    wx[ilayer] += hit_x->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
	    wy[ilayer] += hit_y->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
	  }
	}
	  bary_x[ilayer] = wx[ilayer] / sume_w[ilayer];
	  bary_y[ilayer] = wy[ilayer] / sume_w[ilayer];
	}
      if( (sume_w[ilayer] == 0) or (sume_w[ilayer] < 0) ){
	bary_x[ilayer] = 0.;
	bary_y[ilayer] = 0.;
      }
      blv[ilayer][0] = bary_x[ilayer];
      blv[ilayer][1] = bary_y[ilayer];
      blv[ilayer][2] = pow(pow(bary_x[ilayer],2) + pow(bary_y[ilayer],2),0.5);
    }
  }
  
  return 0;
}

void fit_bar(TH1 * h, double &mu, double &sig){
  mu = 0.;
  sig = 0.;

  if(h->GetEntries() > 10) {
    h->Fit("gaus", "q");
    mu = h->GetFunction("gaus")->GetParameter("Mean");
    sig = h->GetFunction("gaus")->GetParameter("Sigma");
  }
  //cout<<"Fit bar"<<endl;
  //cout<<"Gaus fit value: "<<mu<<" +- "<<sig<<endl;
  return;
}

void fit_res(TH1 * h, double &mu, double &sig, double &res){
    h->Fit("gaus", "q");
    mu = h->GetFunction("gaus")->GetParameter("Mean");
    sig = h->GetFunction("gaus")->GetParameter("Sigma");
    res = h->GetFunction("gaus")->GetParameter("Sigma") / h->GetFunction("gaus")->GetParameter("Mean");
    return;
}

void fit_gauss(TH1 * h, double &mu, double &mu_error, double &sig){

  mu = 0.;
  mu_error = 0.;
  sig = 0.;

  int h_binmax = h->GetMaximumBin();
  int h_NBINS = h->GetNbinsX();
  //double h_xmaxvalue = h->GetXaxis()->GetBinCenter(h_binmax);                                                                                               
  double h_MEAN = h->GetMean();
  double h_maxvalue = h->GetBinContent(h_binmax);
  double h_mpv = h->GetXaxis()->GetBinCenter(h_binmax);
  double h_std = h->GetStdDev();
  double xmin = h->GetXaxis()->GetBinCenter(h_binmax)-2*h_std;
  if(xmin <0.) xmin=0.001;
  double xmax = h->GetXaxis()->GetBinCenter(h_binmax)+2*h_std;

  TF1 *GaussDist = new TF1("GaussDist","[0]*TMath::Exp(-0.5*pow((x-[1])/[2],2))",xmin,xmax);
  GaussDist->SetParameters(0.2*abs(xmax-xmin)*h_maxvalue,h_mpv,h_std);

  h->Fit(GaussDist,"LR");
  mu = GaussDist->GetParameter(1);
  mu_error = GaussDist->GetParError(1);
  sig = GaussDist->GetParameter(2);

  return;

}

double langaufun(double *x, double *par) {
 
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
 
  // Numeric constants
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double mpshift  = -0.22278298;       // Landau maximum location
 
  // Control constants
  double np = 350.0;      // number of convolution steps
  double sc =   3.5;      // convolution extends to +-sc Gaussian sigmas
 
  // Variables
  double xx;
  double mpc;
  double fland;
  double sum = 0.0;
  double xlow,xupp;
  double step;
  double i;
 
 
  // MP shift correction
  mpc = par[1] - mpshift * par[0];
 
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
 
  step = (xupp-xlow) / np;
 
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
 
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
 
  return (par[2] * step * sum * invsq2pi / par[3]);
}

void fit_langaus(TH1 * h, double &MPV, double &MPV_error, double &FWHM){

  MPV = 0.;
  MPV_error = 0.;
  FWHM = 0.;

  int h_binmax = h->GetMaximumBin();
  int h_NBINS = h->GetNbinsX();
  //double h_xmaxvalue = h->GetXaxis()->GetBinCenter(h_binmax);
  double h_MEAN = h->GetMean();
  double h_maxvalue = h->GetBinContent(h_binmax);
  double h_mpv = h->GetXaxis()->GetBinCenter(h_binmax);
  double h_std = h->GetStdDev();
  double xmin = h->GetXaxis()->GetBinCenter(h_binmax)-4*h_std;
  if(xmin < 0) xmin = 0.01;
  double xmax = h->GetXaxis()->GetBinCenter(h_binmax)+4*h_std;
  
  double maxvalue = 0.;
    
  cout<<"FITTING: "<<h->GetName()<<endl;
    if(h_mpv > 0.) { 
      if(h_std > 0.15*h_MEAN){
	cout<<"Landau-Gauss"<<endl;
	TF1 *LandauGaussDist = new TF1("LandauGaussDist",langaufun,xmin,xmax,4);
	LandauGaussDist->SetParameters(h_std,h_mpv,0.2*abs(xmax-xmin)*h_maxvalue,h_std);
	
	h->Fit(LandauGaussDist,"LR");
	cout<<"Parameters: "<<LandauGaussDist->GetParameter(0)<<", "<<LandauGaussDist->GetParameter(1)<<", "<<LandauGaussDist->GetParameter(2)<<", "<<LandauGaussDist->GetParameter(3)<<endl;
    
	MPV = LandauGaussDist->GetMaximumX();
	MPV_error = LandauGaussDist->GetParError(1);
	
	maxvalue = LandauGaussDist->Eval(MPV);
     
	//cout<<"Max. Value: "<<maxvalue<<", at x= "<<MPV<<endl;
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
        
	if((FWHM > 2.5*h_std) || (FWHM < h_std) || (MPV > 1.5*h_MEAN) || (MPV < 0.5*h_MEAN)) {
	  MPV_error = 0.;
	  MPV = h_mpv;
	  FWHM = 2*h_std;
	}
      }
      else{
	// Gaussian
	cout<<"Bypass to Gauss"<<endl;
	TF1 *GaussDist = new TF1("GaussDist","[0]*TMath::Exp(-0.5*pow((x-[1])/[2],2))",xmin,xmax);
        GaussDist->SetParameters(0.2*abs(xmax-xmin)*h_maxvalue,h_mpv,h_std);
		  
	h->Fit(GaussDist,"LR");
	MPV = GaussDist->GetParameter(1);
	MPV_error = GaussDist->GetParError(1);
	FWHM = 2*GaussDist->GetParameter(2);
	cout<<"Parameters: "<<GaussDist->GetParameter(1)<<", "<<GaussDist->GetParameter(2)<<endl;
      }
    }
    else{
      MPV_error = 0.;
      MPV = 0.;
      FWHM = 0.;
    }
    
  cout<<"My value: "<<MPV<<" +- "<<FWHM/2<<endl;

  return; 
}


void graph_setup_add(TGraph *g, string title, Color_t color){
    g->SetTitle(title.c_str());
    g->SetLineColor(color);
    g->SetLineWidth(3);
    return;
}


void analysis_125 (string particle, bool masked=false) {
    
    // double energies[N_ENERGIES] = {1, 2, 5, 8, 10, 20, 40, 60, 80, 100, 150};
    double test_e[N_ENERGIES]={125};
    //double test_e[N_ENERGIES]={2., 150.};
    //double test_e[N_ENERGIES]={2., 40., 150.};
    //double test_e[N_ENERGIES] = {2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200};
    TVectorD energies(N_ENERGIES, test_e);
    TVectorD energies_tr(N_ENERGIES);
    for (int j = 0; j < N_ENERGIES; j++) energies_tr[j] = 1/TMath::Sqrt(energies[j]);
    // For conf11
    double W[N_ECAL_LAYERS] = {4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6}; 
    TVectorD W_thicknesses(N_ECAL_LAYERS, W);
    //Si: 650 650 650 650 500 500 500 500 500 500 320 320 320 320 320
    
    string filenames[N_ENERGIES];
    string base_path = "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/submit_jobs/LCIO2build_folder/LCIO2build_output/";
    // for (int j = 0; j < N_ENERGIES; j++) filenames[j] = base_path + "/CONF11/build/ECAL_QGSP_BERT_conf6_e-_" + to_string((int)round(energies[j])) +  "GeV_5kevt_build_masked.root";
    for (int j = 0; j < N_ENERGIES; j++) filenames[j] = base_path + "output_LCIO2Build_TB2022-06_"+particle+"_" + to_string((int)round(energies[j])) +  "GeV.root";
    
    // for (int j = 0; j < N_ENERGIES; j++) filenames[j] = base_path + "CONF11/build/ECAL_QGSP_BERT_conf8_e-_" + to_string((int)round(energies[j])) +  "GeV_5kevt_build_masked.root";
    //double test_zeros[N_ENERGIES]={0.};
    //double test_zeros[N_ENERGIES] = {0., 0.};
    //double test_zeros[N_ENERGIES]={0., 0., 0.};
    double test_zeros[N_ENERGIES];
    for(int i = 0 ; i < N_ENERGIES ; i++) test_zeros[i] = 0;
    
    TVectorD zeros(N_ENERGIES, test_zeros);
    
    TVectorD mu_nhit(N_ENERGIES), mu_sume(N_ENERGIES), mu_weight(N_ENERGIES);
    TVectorD sig_nhit(N_ENERGIES), sig_sume(N_ENERGIES), sig_weight(N_ENERGIES);
    TVectorD res_nhit(N_ENERGIES), res_sume(N_ENERGIES), res_weight(N_ENERGIES);

    TVectorD mu_mol(N_ENERGIES);
    TVectorD mu_error_mol(N_ENERGIES);
    TVectorD sig_mol(N_ENERGIES);
    

    TVectorD mu_barycenter_x(N_ENERGIES), mu_barycenter_y(N_ENERGIES), mu_barycenter_z(N_ENERGIES);
    TVectorD sig_barycenter_x(N_ENERGIES), sig_barycenter_y(N_ENERGIES), sig_barycenter_z(N_ENERGIES);
    TVectorD mu_barycenter_r(N_ENERGIES), sig_barycenter_r(N_ENERGIES);

    TVectorD mu_barycenter_r_layer_0(N_ENERGIES), sig_barycenter_r_layer_0(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_1(N_ENERGIES), sig_barycenter_r_layer_1(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_2(N_ENERGIES), sig_barycenter_r_layer_2(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_3(N_ENERGIES), sig_barycenter_r_layer_3(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_4(N_ENERGIES), sig_barycenter_r_layer_4(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_5(N_ENERGIES), sig_barycenter_r_layer_5(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_6(N_ENERGIES), sig_barycenter_r_layer_6(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_7(N_ENERGIES), sig_barycenter_r_layer_7(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_8(N_ENERGIES), sig_barycenter_r_layer_8(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_9(N_ENERGIES), sig_barycenter_r_layer_9(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_10(N_ENERGIES), sig_barycenter_r_layer_10(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_11(N_ENERGIES), sig_barycenter_r_layer_11(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_12(N_ENERGIES), sig_barycenter_r_layer_12(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_13(N_ENERGIES), sig_barycenter_r_layer_13(N_ENERGIES);
    TVectorD mu_barycenter_r_layer_14(N_ENERGIES), sig_barycenter_r_layer_14(N_ENERGIES);

    TVectorD mu_barycenter_x_layer_0(N_ENERGIES), sig_barycenter_x_layer_0(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_1(N_ENERGIES), sig_barycenter_x_layer_1(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_2(N_ENERGIES), sig_barycenter_x_layer_2(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_3(N_ENERGIES), sig_barycenter_x_layer_3(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_4(N_ENERGIES), sig_barycenter_x_layer_4(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_5(N_ENERGIES), sig_barycenter_x_layer_5(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_6(N_ENERGIES), sig_barycenter_x_layer_6(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_7(N_ENERGIES), sig_barycenter_x_layer_7(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_8(N_ENERGIES), sig_barycenter_x_layer_8(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_9(N_ENERGIES), sig_barycenter_x_layer_9(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_10(N_ENERGIES), sig_barycenter_x_layer_10(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_11(N_ENERGIES), sig_barycenter_x_layer_11(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_12(N_ENERGIES), sig_barycenter_x_layer_12(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_13(N_ENERGIES), sig_barycenter_x_layer_13(N_ENERGIES);
    TVectorD mu_barycenter_x_layer_14(N_ENERGIES), sig_barycenter_x_layer_14(N_ENERGIES);

    TVectorD mu_barycenter_y_layer_0(N_ENERGIES), sig_barycenter_y_layer_0(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_1(N_ENERGIES), sig_barycenter_y_layer_1(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_2(N_ENERGIES), sig_barycenter_y_layer_2(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_3(N_ENERGIES), sig_barycenter_y_layer_3(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_4(N_ENERGIES), sig_barycenter_y_layer_4(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_5(N_ENERGIES), sig_barycenter_y_layer_5(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_6(N_ENERGIES), sig_barycenter_y_layer_6(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_7(N_ENERGIES), sig_barycenter_y_layer_7(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_8(N_ENERGIES), sig_barycenter_y_layer_8(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_9(N_ENERGIES), sig_barycenter_y_layer_9(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_10(N_ENERGIES), sig_barycenter_y_layer_10(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_11(N_ENERGIES), sig_barycenter_y_layer_11(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_12(N_ENERGIES), sig_barycenter_y_layer_12(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_13(N_ENERGIES), sig_barycenter_y_layer_13(N_ENERGIES);
    TVectorD mu_barycenter_y_layer_14(N_ENERGIES), sig_barycenter_y_layer_14(N_ENERGIES);

    TVectorD mu_nhit_layer_0(N_ENERGIES), mu_error_nhit_layer_0(N_ENERGIES), sig_nhit_layer_0(N_ENERGIES);
    TVectorD mu_nhit_layer_1(N_ENERGIES), mu_error_nhit_layer_1(N_ENERGIES), sig_nhit_layer_1(N_ENERGIES);
    TVectorD mu_nhit_layer_2(N_ENERGIES), mu_error_nhit_layer_2(N_ENERGIES), sig_nhit_layer_2(N_ENERGIES);
    TVectorD mu_nhit_layer_3(N_ENERGIES), mu_error_nhit_layer_3(N_ENERGIES), sig_nhit_layer_3(N_ENERGIES);
    TVectorD mu_nhit_layer_4(N_ENERGIES), mu_error_nhit_layer_4(N_ENERGIES), sig_nhit_layer_4(N_ENERGIES);
    TVectorD mu_nhit_layer_5(N_ENERGIES), mu_error_nhit_layer_5(N_ENERGIES), sig_nhit_layer_5(N_ENERGIES);
    TVectorD mu_nhit_layer_6(N_ENERGIES), mu_error_nhit_layer_6(N_ENERGIES), sig_nhit_layer_6(N_ENERGIES);
    TVectorD mu_nhit_layer_7(N_ENERGIES), mu_error_nhit_layer_7(N_ENERGIES), sig_nhit_layer_7(N_ENERGIES);
    TVectorD mu_nhit_layer_8(N_ENERGIES), mu_error_nhit_layer_8(N_ENERGIES), sig_nhit_layer_8(N_ENERGIES);
    TVectorD mu_nhit_layer_9(N_ENERGIES), mu_error_nhit_layer_9(N_ENERGIES), sig_nhit_layer_9(N_ENERGIES);
    TVectorD mu_nhit_layer_10(N_ENERGIES), mu_error_nhit_layer_10(N_ENERGIES), sig_nhit_layer_10(N_ENERGIES);
    TVectorD mu_nhit_layer_11(N_ENERGIES), mu_error_nhit_layer_11(N_ENERGIES), sig_nhit_layer_11(N_ENERGIES);
    TVectorD mu_nhit_layer_12(N_ENERGIES), mu_error_nhit_layer_12(N_ENERGIES), sig_nhit_layer_12(N_ENERGIES);
    TVectorD mu_nhit_layer_13(N_ENERGIES), mu_error_nhit_layer_13(N_ENERGIES), sig_nhit_layer_13(N_ENERGIES);
    TVectorD mu_nhit_layer_14(N_ENERGIES), mu_error_nhit_layer_14(N_ENERGIES), sig_nhit_layer_14(N_ENERGIES);

    TVectorD mu_nhit_layer_n_0(N_ENERGIES), mu_error_nhit_layer_n_0(N_ENERGIES), sig_nhit_layer_n_0(N_ENERGIES);
    TVectorD mu_nhit_layer_n_1(N_ENERGIES), mu_error_nhit_layer_n_1(N_ENERGIES), sig_nhit_layer_n_1(N_ENERGIES);
    TVectorD mu_nhit_layer_n_2(N_ENERGIES), mu_error_nhit_layer_n_2(N_ENERGIES), sig_nhit_layer_n_2(N_ENERGIES);
    TVectorD mu_nhit_layer_n_3(N_ENERGIES), mu_error_nhit_layer_n_3(N_ENERGIES), sig_nhit_layer_n_3(N_ENERGIES);
    TVectorD mu_nhit_layer_n_4(N_ENERGIES), mu_error_nhit_layer_n_4(N_ENERGIES), sig_nhit_layer_n_4(N_ENERGIES);
    TVectorD mu_nhit_layer_n_5(N_ENERGIES), mu_error_nhit_layer_n_5(N_ENERGIES), sig_nhit_layer_n_5(N_ENERGIES);
    TVectorD mu_nhit_layer_n_6(N_ENERGIES), mu_error_nhit_layer_n_6(N_ENERGIES), sig_nhit_layer_n_6(N_ENERGIES);
    TVectorD mu_nhit_layer_n_7(N_ENERGIES), mu_error_nhit_layer_n_7(N_ENERGIES), sig_nhit_layer_n_7(N_ENERGIES);
    TVectorD mu_nhit_layer_n_8(N_ENERGIES), mu_error_nhit_layer_n_8(N_ENERGIES), sig_nhit_layer_n_8(N_ENERGIES);
    TVectorD mu_nhit_layer_n_9(N_ENERGIES), mu_error_nhit_layer_n_9(N_ENERGIES), sig_nhit_layer_n_9(N_ENERGIES);
    TVectorD mu_nhit_layer_n_10(N_ENERGIES), mu_error_nhit_layer_n_10(N_ENERGIES), sig_nhit_layer_n_10(N_ENERGIES);
    TVectorD mu_nhit_layer_n_11(N_ENERGIES), mu_error_nhit_layer_n_11(N_ENERGIES), sig_nhit_layer_n_11(N_ENERGIES);
    TVectorD mu_nhit_layer_n_12(N_ENERGIES), mu_error_nhit_layer_n_12(N_ENERGIES), sig_nhit_layer_n_12(N_ENERGIES);
    TVectorD mu_nhit_layer_n_13(N_ENERGIES), mu_error_nhit_layer_n_13(N_ENERGIES), sig_nhit_layer_n_13(N_ENERGIES);
    TVectorD mu_nhit_layer_n_14(N_ENERGIES), mu_error_nhit_layer_n_14(N_ENERGIES), sig_nhit_layer_n_14(N_ENERGIES);

    TVectorD mu_weight_layer_0(N_ENERGIES), mu_error_weight_layer_0(N_ENERGIES), sig_weight_layer_0(N_ENERGIES);
    TVectorD mu_weight_layer_1(N_ENERGIES), mu_error_weight_layer_1(N_ENERGIES), sig_weight_layer_1(N_ENERGIES);
    TVectorD mu_weight_layer_2(N_ENERGIES), mu_error_weight_layer_2(N_ENERGIES), sig_weight_layer_2(N_ENERGIES);
    TVectorD mu_weight_layer_3(N_ENERGIES), mu_error_weight_layer_3(N_ENERGIES), sig_weight_layer_3(N_ENERGIES);
    TVectorD mu_weight_layer_4(N_ENERGIES), mu_error_weight_layer_4(N_ENERGIES), sig_weight_layer_4(N_ENERGIES);
    TVectorD mu_weight_layer_5(N_ENERGIES), mu_error_weight_layer_5(N_ENERGIES), sig_weight_layer_5(N_ENERGIES);
    TVectorD mu_weight_layer_6(N_ENERGIES), mu_error_weight_layer_6(N_ENERGIES), sig_weight_layer_6(N_ENERGIES);
    TVectorD mu_weight_layer_7(N_ENERGIES), mu_error_weight_layer_7(N_ENERGIES), sig_weight_layer_7(N_ENERGIES);
    TVectorD mu_weight_layer_8(N_ENERGIES), mu_error_weight_layer_8(N_ENERGIES), sig_weight_layer_8(N_ENERGIES);
    TVectorD mu_weight_layer_9(N_ENERGIES), mu_error_weight_layer_9(N_ENERGIES), sig_weight_layer_9(N_ENERGIES);
    TVectorD mu_weight_layer_10(N_ENERGIES), mu_error_weight_layer_10(N_ENERGIES), sig_weight_layer_10(N_ENERGIES);
    TVectorD mu_weight_layer_11(N_ENERGIES), mu_error_weight_layer_11(N_ENERGIES), sig_weight_layer_11(N_ENERGIES);
    TVectorD mu_weight_layer_12(N_ENERGIES), mu_error_weight_layer_12(N_ENERGIES), sig_weight_layer_12(N_ENERGIES);
    TVectorD mu_weight_layer_13(N_ENERGIES), mu_error_weight_layer_13(N_ENERGIES), sig_weight_layer_13(N_ENERGIES);
    TVectorD mu_weight_layer_14(N_ENERGIES), mu_error_weight_layer_14(N_ENERGIES), sig_weight_layer_14(N_ENERGIES);

    TVectorD mu_weight_layer_n_0(N_ENERGIES), mu_error_weight_layer_n_0(N_ENERGIES), sig_weight_layer_n_0(N_ENERGIES);
    TVectorD mu_weight_layer_n_1(N_ENERGIES), mu_error_weight_layer_n_1(N_ENERGIES), sig_weight_layer_n_1(N_ENERGIES);
    TVectorD mu_weight_layer_n_2(N_ENERGIES), mu_error_weight_layer_n_2(N_ENERGIES), sig_weight_layer_n_2(N_ENERGIES);
    TVectorD mu_weight_layer_n_3(N_ENERGIES), mu_error_weight_layer_n_3(N_ENERGIES), sig_weight_layer_n_3(N_ENERGIES);
    TVectorD mu_weight_layer_n_4(N_ENERGIES), mu_error_weight_layer_n_4(N_ENERGIES), sig_weight_layer_n_4(N_ENERGIES);
    TVectorD mu_weight_layer_n_5(N_ENERGIES), mu_error_weight_layer_n_5(N_ENERGIES), sig_weight_layer_n_5(N_ENERGIES);
    TVectorD mu_weight_layer_n_6(N_ENERGIES), mu_error_weight_layer_n_6(N_ENERGIES), sig_weight_layer_n_6(N_ENERGIES);
    TVectorD mu_weight_layer_n_7(N_ENERGIES), mu_error_weight_layer_n_7(N_ENERGIES), sig_weight_layer_n_7(N_ENERGIES);
    TVectorD mu_weight_layer_n_8(N_ENERGIES), mu_error_weight_layer_n_8(N_ENERGIES), sig_weight_layer_n_8(N_ENERGIES);
    TVectorD mu_weight_layer_n_9(N_ENERGIES), mu_error_weight_layer_n_9(N_ENERGIES), sig_weight_layer_n_9(N_ENERGIES);
    TVectorD mu_weight_layer_n_10(N_ENERGIES), mu_error_weight_layer_n_10(N_ENERGIES), sig_weight_layer_n_10(N_ENERGIES);
    TVectorD mu_weight_layer_n_11(N_ENERGIES), mu_error_weight_layer_n_11(N_ENERGIES), sig_weight_layer_n_11(N_ENERGIES);
    TVectorD mu_weight_layer_n_12(N_ENERGIES), mu_error_weight_layer_n_12(N_ENERGIES), sig_weight_layer_n_12(N_ENERGIES);
    TVectorD mu_weight_layer_n_13(N_ENERGIES), mu_error_weight_layer_n_13(N_ENERGIES), sig_weight_layer_n_13(N_ENERGIES);
    TVectorD mu_weight_layer_n_14(N_ENERGIES), mu_error_weight_layer_n_14(N_ENERGIES), sig_weight_layer_n_14(N_ENERGIES);

    TVectorD shower_nhit_max_layer(N_ENERGIES), shower_nhit_start_layer(N_ENERGIES), shower_nhit_end_layer(N_ENERGIES), shower_nhit_start_10_layer(N_ENERGIES), shower_nhit_end_10_layer(N_ENERGIES), shower_nhit_average(N_ENERGIES), shower_nhit_max(N_ENERGIES);
    TVectorD shower_weight_max_layer(N_ENERGIES), shower_weight_start_layer(N_ENERGIES), shower_weight_end_layer(N_ENERGIES), shower_weight_start_10_layer(N_ENERGIES), shower_weight_end_10_layer(N_ENERGIES), shower_weight_average(N_ENERGIES), shower_weight_max(N_ENERGIES);

    // Output filename
    TString result_name = "resolution_"+particle+"_125_result.root" ;
    TFile f(result_name, "recreate");

    // Begin of reading files and writing histos
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
	
        TH1F *h_nhit = new TH1F(("NumHits_" + part_string + "_" + e_str).c_str(), ("Number of hits " + part_string + " " + e_str).c_str(), 4001, -0.5, 4000.5);
	TH1F *h_sume = new TH1F(("SumEnergy_" + part_string + "_" + e_str).c_str(), ("Sum Energy " + part_string + " " + e_str).c_str(), 5001, -0.5, 15000.5);
	TH1F *h_weight = new TH1F(("WSumEnergy_" + part_string + "_" + e_str).c_str(), ("W Sum Energy " + part_string + " " + e_str).c_str(), 5001, -0.5, 50000.5);
                
	// Barycenter histos
        TH1F *h_bar_x = new TH1F(("Barycenter_x_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_y = new TH1F(("Barycenter_y_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_z = new TH1F(("Barycenter_z_" + part_string + "_" + e_str).c_str(), ("Barycenter (z-axis) " + part_string + " " + e_str).c_str(), 1000, -500., 500.);
	TH1F *h_bar_r = new TH1F(("Barycenter_r_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);

	// Barycenter per layer
	TH1F *h_bar_r_layer_0 = new TH1F(("Barycenter_r_layer_0_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 0 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_1 = new TH1F(("Barycenter_r_layer_1_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 1 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_2 = new TH1F(("Barycenter_r_layer_2_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 2 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_3 = new TH1F(("Barycenter_r_layer_3_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 3 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_r_layer_4 = new TH1F(("Barycenter_r_layer_4_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 4 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_5 = new TH1F(("Barycenter_r_layer_5_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 5 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_6 = new TH1F(("Barycenter_r_layer_6_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 6 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_7 = new TH1F(("Barycenter_r_layer_7_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 7 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_8 = new TH1F(("Barycenter_r_layer_8_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 8 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_9 = new TH1F(("Barycenter_r_layer_9_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 9 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_10 = new TH1F(("Barycenter_r_layer_10_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 10 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_11 = new TH1F(("Barycenter_r_layer_11_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 11 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_12 = new TH1F(("Barycenter_r_layer_12_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 12 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_13 = new TH1F(("Barycenter_r_layer_13_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 13 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_r_layer_14 = new TH1F(("Barycenter_r_layer_14_" + part_string + "_" + e_str).c_str(), ("Barycenter (radial) layer 14 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);

	TH1F *h_bar_x_layer_0 = new TH1F(("Barycenter_x_layer_0_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 0 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_1 = new TH1F(("Barycenter_x_layer_1_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 1 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_2 = new TH1F(("Barycenter_x_layer_2_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 2 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_3 = new TH1F(("Barycenter_x_layer_3_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 3 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_4 = new TH1F(("Barycenter_x_layer_4_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 4 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_5 = new TH1F(("Barycenter_x_layer_5_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 5 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_6 = new TH1F(("Barycenter_x_layer_6_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 6 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_7 = new TH1F(("Barycenter_x_layer_7_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 7 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_8 = new TH1F(("Barycenter_x_layer_8_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 8 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_9 = new TH1F(("Barycenter_x_layer_9_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 9 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_10 = new TH1F(("Barycenter_x_layer_10_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 10 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_11 = new TH1F(("Barycenter_x_layer_11_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 11 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_12 = new TH1F(("Barycenter_x_layer_12_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 12 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_13 = new TH1F(("Barycenter_x_layer_13_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 13 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_x_layer_14 = new TH1F(("Barycenter_x_layer_14_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 14 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);

	TH1F *h_bar_y_layer_0 = new TH1F(("Barycenter_y_layer_0_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 0 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_y_layer_1 = new TH1F(("Barycenter_y_layer_1_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 1 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_y_layer_2 = new TH1F(("Barycenter_y_layer_2_" + part_string + "_" + e_str).c_str(), ("Barycenter (x-axis) layer 2 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_y_layer_3 = new TH1F(("Barycenter_y_layer_3_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 3 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_y_layer_4 = new TH1F(("Barycenter_y_layer_4_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 4 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_y_layer_5 = new TH1F(("Barycenter_y_layer_5_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 5 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_y_layer_6 = new TH1F(("Barycenter_y_layer_6_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 6 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_y_layer_7 = new TH1F(("Barycenter_y_layer_7_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 7 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_y_layer_8 = new TH1F(("Barycenter_y_layer_8_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 8 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
        TH1F *h_bar_y_layer_9 = new TH1F(("Barycenter_y_layer_9_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 9 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_y_layer_10 = new TH1F(("Barycenter_y_layer_10_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 10 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_y_layer_11 = new TH1F(("Barycenter_y_layer_11_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 11 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_y_layer_12 = new TH1F(("Barycenter_y_layer_12_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 12 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_y_layer_13 = new TH1F(("Barycenter_y_layer_13_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 13 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);
	TH1F *h_bar_y_layer_14 = new TH1F(("Barycenter_y_layer_14_" + part_string + "_" + e_str).c_str(), ("Barycenter (y-axis) layer 14 " + part_string + " " + e_str).c_str(), 201, -100.5, 100.5);

	// Moliere histos
        TH1F *h_mol = new TH1F(("Radius90_" + part_string + "_" + e_str).c_str(), ("Radius containing 90% of energy " + part_string + " " + e_str).c_str(), 101, -0.5, 100.5);
        
	// Shower profile characteristics
	TH1F *h_shower_nhit_max_layer = new TH1F(("ShowerNhitMaxLayer_" + part_string + "_" + e_str).c_str(), ("Shower Nhit Max. (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
	TH1F *h_shower_nhit_start_layer = new TH1F(("ShowerNhitStartLayer_" + part_string + "_" + e_str).c_str(), ("Shower Nhit Start (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
	TH1F *h_shower_nhit_end_layer = new TH1F(("ShowerNhitEndLayer_" + part_string + "_" + e_str).c_str(), ("Shower Nhit End (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
	TH1F *h_shower_nhit_start_10_layer = new TH1F(("ShowerNhitStart10Layer_" + part_string + "_" + e_str).c_str(), ("Shower Nhit Start 0.1*Max (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
        TH1F *h_shower_nhit_end_10_layer = new TH1F(("ShowerNhitEnd10Layer_" + part_string + "_" + e_str).c_str(), ("Shower Nhit End 0.1*Max (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
	TH1F *h_shower_nhit_average= new TH1F(("ShowerNhitAverage_" + part_string + "_" + e_str).c_str(), ("Shower Nhit Average  " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_shower_nhit_max= new TH1F(("ShowerNhitMax_" + part_string + "_" + e_str).c_str(), ("Shower Nhit Max  " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);

	TH1F *h_shower_weight_max_layer = new TH1F(("ShowerWeightMaxLayer_" + part_string + "_" + e_str).c_str(), ("Shower Weight Max. (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
        TH1F *h_shower_weight_start_layer = new TH1F(("ShowerWeightStartLayer_" + part_string + "_" + e_str).c_str(), ("Shower Weight Start (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
        TH1F *h_shower_weight_end_layer = new TH1F(("ShowerWeightEndLayer_" + part_string + "_" + e_str).c_str(), ("Shower Weight End (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
        TH1F *h_shower_weight_start_10_layer = new TH1F(("ShowerWeightStart10Layer_" + part_string + "_" + e_str).c_str(), ("Shower Weight Start 0.1*Max (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
        TH1F *h_shower_weight_end_10_layer = new TH1F(("ShowerWeightEnd10Layer_" + part_string + "_" + e_str).c_str(), ("Shower Weight End 0.1*Max (layer) " + part_string + " " + e_str).c_str(), 16, -1.5, 14.5);
	TH1F *h_shower_weight_average= new TH1F(("ShowerWeightAverage_" + part_string + "_" + e_str).c_str(), ("Shower Weight Average  " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_shower_weight_max= new TH1F(("ShowerWeightMax_" + part_string + "_" + e_str).c_str(), ("Shower Weight Max  " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
	
	// Shower profile per layer
	TH1F *h_nhit_layer_0 = new TH1F(("NumHits_layer_0_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 0) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_1 = new TH1F(("NumHits_layer_1_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 1) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_2 = new TH1F(("NumHits_layer_2_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 2) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_3 = new TH1F(("NumHits_layer_3_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 3) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_4 = new TH1F(("NumHits_layer_4_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 4) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_5 = new TH1F(("NumHits_layer_5_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 5) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_6 = new TH1F(("NumHits_layer_6_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 6) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_7 = new TH1F(("NumHits_layer_7_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 7) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_8 = new TH1F(("NumHits_layer_8_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 8) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_9 = new TH1F(("NumHits_layer_9_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 9) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_10 = new TH1F(("NumHits_layer_10_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 10) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_11 = new TH1F(("NumHits_layer_11_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 11) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_12 = new TH1F(("NumHits_layer_12_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 12) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_13 = new TH1F(("NumHits_layer_13_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 13) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
        TH1F *h_nhit_layer_14 = new TH1F(("NumHits_layer_14_" + part_string + "_" + e_str).c_str(), ("Number of hits (layer 14) " + part_string + " " + e_str).c_str(), 201, -0.5, 200.5);
	
	TH1F *h_nhit_layer_n_0 = new TH1F(("NumHits_layer_n_0_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 0) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_1 = new TH1F(("NumHits_layer_n_1_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 1) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_2 = new TH1F(("NumHits_layer_n_2_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 2) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_3 = new TH1F(("NumHits_layer_n_3_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 3) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_4 = new TH1F(("NumHits_layer_n_4_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 4) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_5 = new TH1F(("NumHits_layer_n_5_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 5) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_6 = new TH1F(("NumHits_layer_n_6_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 6) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_7 = new TH1F(("NumHits_layer_n_7_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 7) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_8 = new TH1F(("NumHits_layer_n_8_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 8) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_9 = new TH1F(("NumHits_layer_n_9_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 9) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_10 = new TH1F(("NumHits_layer_n_10_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 10) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_11 = new TH1F(("NumHits_layer_n_11_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 11) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_12 = new TH1F(("NumHits_layer_n_12_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 12) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_13 = new TH1F(("NumHits_layer_n_13_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 13) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_nhit_layer_n_14 = new TH1F(("NumHits_layer_n_14_" + part_string + "_" + e_str).c_str(), ("Number of hits normalized (layer 14) " + part_string + " " + e_str).c_str(), 100, 0, 1);
		
	TH1F *h_weight_layer_0 = new TH1F(("Weight_layer_0_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 0) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_1 = new TH1F(("Weight_layer_1_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 1) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_2 = new TH1F(("Weight_layer_2_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 2) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_3 = new TH1F(("Weight_layer_3_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 3) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_4 = new TH1F(("Weight_layer_4_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 4) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_5 = new TH1F(("Weight_layer_5_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 5) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_6 = new TH1F(("Weight_layer_6_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 6) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_7 = new TH1F(("Weight_layer_7_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 7) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_8 = new TH1F(("Weight_layer_8_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 8) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_9 = new TH1F(("Weight_layer_9_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 9) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_10 = new TH1F(("Weight_layer_10_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 10) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_11 = new TH1F(("Weight_layer_11_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 11) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_12 = new TH1F(("Weight_layer_12_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 12) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_13 = new TH1F(("Weight_layer_13_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 13) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        TH1F *h_weight_layer_14 = new TH1F(("Weight_layer_14_" + part_string + "_" + e_str).c_str(), ("Weighted energy (layer 14) " + part_string + " " + e_str).c_str(), 1001, -0.5, 10000.5);
        
	TH1F *h_weight_layer_n_0 = new TH1F(("Weight_layer_n_0_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 0) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_1 = new TH1F(("Weight_layer_n_1_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 1) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_2 = new TH1F(("Weight_layer_n_2_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 2) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_3 = new TH1F(("Weight_layer_n_3_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 3) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_4 = new TH1F(("Weight_layer_n_4_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 4) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_5 = new TH1F(("Weight_layer_n_5_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 5) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_6 = new TH1F(("Weight_layer_n_6_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 6) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_7 = new TH1F(("Weight_layer_n_7_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 7) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_8 = new TH1F(("Weight_layer_n_8_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 8) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_9 = new TH1F(("Weight_layer_n_9_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 9) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_10 = new TH1F(("Weight_layer_n_10_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 10) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_11 = new TH1F(("Weight_layer_n_11_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 11) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_12 = new TH1F(("Weight_layer_n_12_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 12) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_13 = new TH1F(("Weight_layer_n_13_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 13) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        TH1F *h_weight_layer_n_14 = new TH1F(("Weight_layer_n_14_" + part_string + "_" + e_str).c_str(), ("Weighted energy normalized (layer 14) " + part_string + " " + e_str).c_str(), 100, 0, 1);
        	
	for (int i_event = 0; i_event < (int)nentries; i_event++) {
	  // Resolution
	  int nhit = 0;     
	  float sume = 0;   
	  float weight = 0; 
	  
	  float bar_xyzr[4]; // 4 for (x,y,z,r)
	  
	  float nhit_layer_array[N_ECAL_LAYERS];
	  float nhit_layer_n_array[N_ECAL_LAYERS];

	  float weight_layer_array[N_ECAL_LAYERS];
          float weight_layer_n_array[N_ECAL_LAYERS];

	  float bar_layer_array[N_ECAL_LAYERS][3]; // 3 for (x,y,z)

	  tree->GetEntry(i_event);
	  if(i_event % 1000 == 0) cout << "Event " << to_string(i_event) << endl;
	  
	  get_res(nhit,	sume, weight, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked);

	  h_nhit->Fill(nhit);        
	  h_sume->Fill(sume);         
	  h_weight->Fill(weight);     
	  
	  // Fill barycenter
	  barycenter(hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, bar_xyzr, hit_isMasked, masked);
	  
	  h_bar_x->Fill(bar_xyzr[0]);
	  h_bar_y->Fill(bar_xyzr[1]);
	  h_bar_z->Fill(bar_xyzr[2]);
	  h_bar_r->Fill(bar_xyzr[3]);
	  
	  // Fill Moliere radii histograms
	  h_mol->Fill(moliere(hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, hit_isMasked, masked));
	  
	  // Fill shower profile 
	  // Nhit
	  hits_layer(nhit_layer_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, false, "nhit");
	  hits_layer(nhit_layer_n_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, true, "nhit");
          
	  h_nhit_layer_0->Fill(nhit_layer_array[0]);
	  h_nhit_layer_1->Fill(nhit_layer_array[1]);
	  h_nhit_layer_2->Fill(nhit_layer_array[2]);
	  h_nhit_layer_3->Fill(nhit_layer_array[3]);
	  h_nhit_layer_4->Fill(nhit_layer_array[4]);
	  h_nhit_layer_5->Fill(nhit_layer_array[5]);
	  h_nhit_layer_6->Fill(nhit_layer_array[6]);
	  h_nhit_layer_7->Fill(nhit_layer_array[7]);
	  h_nhit_layer_8->Fill(nhit_layer_array[8]);
	  h_nhit_layer_9->Fill(nhit_layer_array[9]);
	  h_nhit_layer_10->Fill(nhit_layer_array[10]);
	  h_nhit_layer_11->Fill(nhit_layer_array[11]);
	  h_nhit_layer_12->Fill(nhit_layer_array[12]);
	  h_nhit_layer_13->Fill(nhit_layer_array[13]);
	  h_nhit_layer_14->Fill(nhit_layer_array[14]);
	  	  
	  h_nhit_layer_n_0->Fill(nhit_layer_n_array[0]);
	  h_nhit_layer_n_1->Fill(nhit_layer_n_array[1]);
	  h_nhit_layer_n_2->Fill(nhit_layer_n_array[2]);
	  h_nhit_layer_n_3->Fill(nhit_layer_n_array[3]);
	  h_nhit_layer_n_4->Fill(nhit_layer_n_array[4]);
	  h_nhit_layer_n_5->Fill(nhit_layer_n_array[5]);
	  h_nhit_layer_n_6->Fill(nhit_layer_n_array[6]);
	  h_nhit_layer_n_7->Fill(nhit_layer_n_array[7]);
	  h_nhit_layer_n_8->Fill(nhit_layer_n_array[8]);
	  h_nhit_layer_n_9->Fill(nhit_layer_n_array[9]);
	  h_nhit_layer_n_10->Fill(nhit_layer_n_array[10]);
	  h_nhit_layer_n_11->Fill(nhit_layer_n_array[11]);
	  h_nhit_layer_n_12->Fill(nhit_layer_n_array[12]);
	  h_nhit_layer_n_13->Fill(nhit_layer_n_array[13]);
	  h_nhit_layer_n_14->Fill(nhit_layer_n_array[14]);
	  
	  	  
	  // shower nhit general parameters start
	  float nhit_shower_maxvalue = 0.;
	  float nhit_shower_maxvalue_n = 0.;
	  int nhit_ilayermax = -1;
	  int nhit_ilayerstart = -1;
	  int nhit_ilayerend = -1;
	  int nhit_ilayerstart_10 = -1;
          int nhit_ilayerend_10 = -1;
	  
	  float nhit_shower_averagevalue = nhit/N_ECAL_LAYERS;

	  shower_variables(nhit, nhit_layer_array, nhit_layer_n_array, nhit_shower_maxvalue, nhit_shower_maxvalue_n, nhit_ilayermax,
			   nhit_ilayerstart, nhit_ilayerstart_10, nhit_ilayerend, nhit_ilayerend_10, "nhit");
	  //cout<<"nhit shower variables: maxvalue, maxvalue_n,ilayermax, ilayerstart, ilayerstart(10%), ilayerend, ilayerend(10%)"<<endl;
	  //cout<<"nhit shower variables: "<<nhit<<" "<<nhit_shower_maxvalue<<" "<<nhit_shower_maxvalue_n<<" "<<nhit_ilayermax<<" "<<
	  //  nhit_ilayerstart<<" "<<nhit_ilayerstart_10<<" "<<nhit_ilayerend<<" "<<nhit_ilayerend_10<<endl;

	  //shower nhit general parameters finish
	  
	  h_shower_nhit_max_layer->Fill(nhit_ilayermax);
	  h_shower_nhit_start_layer->Fill(nhit_ilayerstart);
	  h_shower_nhit_end_layer->Fill(nhit_ilayerend);
	  h_shower_nhit_start_10_layer->Fill(nhit_ilayerstart_10);
          h_shower_nhit_end_10_layer->Fill(nhit_ilayerend_10);
	  h_shower_nhit_average->Fill(nhit_shower_averagevalue);
	  h_shower_nhit_max->Fill(nhit_shower_maxvalue);
	  
	  // Weighted energy
          hits_layer(weight_layer_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, false, "weight");
	  hits_layer(weight_layer_n_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, true, "weight");
          
          h_weight_layer_0->Fill(weight_layer_array[0]);
          h_weight_layer_1->Fill(weight_layer_array[1]);
          h_weight_layer_2->Fill(weight_layer_array[2]);
          h_weight_layer_3->Fill(weight_layer_array[3]);
          h_weight_layer_4->Fill(weight_layer_array[4]);
          h_weight_layer_5->Fill(weight_layer_array[5]);
          h_weight_layer_6->Fill(weight_layer_array[6]);
          h_weight_layer_7->Fill(weight_layer_array[7]);
          h_weight_layer_8->Fill(weight_layer_array[8]);
          h_weight_layer_9->Fill(weight_layer_array[9]);
          h_weight_layer_10->Fill(weight_layer_array[10]);
          h_weight_layer_11->Fill(weight_layer_array[11]);
          h_weight_layer_12->Fill(weight_layer_array[12]);
          h_weight_layer_13->Fill(weight_layer_array[13]);
          h_weight_layer_14->Fill(weight_layer_array[14]);

	  h_weight_layer_n_0->Fill(weight_layer_n_array[0]);
          h_weight_layer_n_1->Fill(weight_layer_n_array[1]);
          h_weight_layer_n_2->Fill(weight_layer_n_array[2]);
          h_weight_layer_n_3->Fill(weight_layer_n_array[3]);
          h_weight_layer_n_4->Fill(weight_layer_n_array[4]);
          h_weight_layer_n_5->Fill(weight_layer_n_array[5]);
          h_weight_layer_n_6->Fill(weight_layer_n_array[6]);
          h_weight_layer_n_7->Fill(weight_layer_n_array[7]);
          h_weight_layer_n_8->Fill(weight_layer_n_array[8]);
          h_weight_layer_n_9->Fill(weight_layer_n_array[9]);
          h_weight_layer_n_10->Fill(weight_layer_n_array[10]);
          h_weight_layer_n_11->Fill(weight_layer_n_array[11]);
          h_weight_layer_n_12->Fill(weight_layer_n_array[12]);
          h_weight_layer_n_13->Fill(weight_layer_n_array[13]);
          h_weight_layer_n_14->Fill(weight_layer_n_array[14]);

	  // shower weight general parameters start
          float weight_shower_maxvalue = 0.;
	  float weight_shower_maxvalue_n = 0.;
          int weight_ilayermax = -1;
          int weight_ilayerstart = -1;
          int weight_ilayerend = -1;
          int weight_ilayerstart_10 = -1;
          int weight_ilayerend_10 = -1;

	  float weight_shower_averagevalue = weight/N_ECAL_LAYERS;

	  shower_variables(weight, weight_layer_array, weight_layer_n_array, weight_shower_maxvalue, weight_shower_maxvalue_n, weight_ilayermax,
                           weight_ilayerstart, weight_ilayerstart_10, weight_ilayerend, weight_ilayerend_10, "weight");
	  //cout<<"weight shower variables: maxvalue, maxvalue_n,ilayermax, ilayerstart, ilayerstart(10%), ilayerend, ilayerend(10%)"<<endl;
	  //cout<<"weight shower variables: "<<weight<<" "<<weight_shower_maxvalue<<" "<<weight_shower_maxvalue_n<<" "<<weight_ilayermax<<" "<<
          //  weight_ilayerstart<<" "<<weight_ilayerstart_10<<" "<<weight_ilayerend<<" "<<weight_ilayerend_10<<endl;

	  //shower weight general parameters finish
	  
          h_shower_weight_max_layer->Fill(weight_ilayermax);
          h_shower_weight_start_layer->Fill(weight_ilayerstart);
          h_shower_weight_end_layer->Fill(weight_ilayerend);
          h_shower_weight_start_10_layer->Fill(weight_ilayerstart_10);
          h_shower_weight_end_10_layer->Fill(weight_ilayerend_10);
          h_shower_weight_average->Fill(weight_shower_averagevalue);
          h_shower_weight_max->Fill(weight_shower_maxvalue);
	  
	  // Barycenter x and y
	  bary_layer(bar_layer_array, hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, hit_isMasked, masked);
	  
          h_bar_x_layer_0->Fill(bar_layer_array[0][0]);
          h_bar_x_layer_1->Fill(bar_layer_array[1][0]);
          h_bar_x_layer_2->Fill(bar_layer_array[2][0]);
          h_bar_x_layer_3->Fill(bar_layer_array[3][0]);
          h_bar_x_layer_4->Fill(bar_layer_array[4][0]);
          h_bar_x_layer_5->Fill(bar_layer_array[5][0]);
          h_bar_x_layer_6->Fill(bar_layer_array[6][0]);
          h_bar_x_layer_7->Fill(bar_layer_array[7][0]);
          h_bar_x_layer_8->Fill(bar_layer_array[8][0]);
          h_bar_x_layer_9->Fill(bar_layer_array[9][0]);
          h_bar_x_layer_10->Fill(bar_layer_array[10][0]);
          h_bar_x_layer_11->Fill(bar_layer_array[11][0]);
          h_bar_x_layer_12->Fill(bar_layer_array[12][0]);
          h_bar_x_layer_13->Fill(bar_layer_array[13][0]);
          h_bar_x_layer_14->Fill(bar_layer_array[14][0]);

	  h_bar_y_layer_0->Fill(bar_layer_array[0][1]);
          h_bar_y_layer_1->Fill(bar_layer_array[1][1]);
          h_bar_y_layer_2->Fill(bar_layer_array[2][1]);
          h_bar_y_layer_3->Fill(bar_layer_array[3][1]);
          h_bar_y_layer_4->Fill(bar_layer_array[4][1]);
          h_bar_y_layer_5->Fill(bar_layer_array[5][1]);
          h_bar_y_layer_6->Fill(bar_layer_array[6][1]);
          h_bar_y_layer_7->Fill(bar_layer_array[7][1]);
          h_bar_y_layer_8->Fill(bar_layer_array[8][1]);
          h_bar_y_layer_9->Fill(bar_layer_array[9][1]);
          h_bar_y_layer_10->Fill(bar_layer_array[10][1]);
          h_bar_y_layer_11->Fill(bar_layer_array[11][1]);
          h_bar_y_layer_12->Fill(bar_layer_array[12][1]);
          h_bar_y_layer_13->Fill(bar_layer_array[13][1]);
          h_bar_y_layer_14->Fill(bar_layer_array[14][1]);

	  h_bar_r_layer_0->Fill(bar_layer_array[0][2]);
          h_bar_r_layer_1->Fill(bar_layer_array[1][2]);
          h_bar_r_layer_2->Fill(bar_layer_array[2][2]);
          h_bar_r_layer_3->Fill(bar_layer_array[3][2]);
          h_bar_r_layer_4->Fill(bar_layer_array[4][2]);
          h_bar_r_layer_5->Fill(bar_layer_array[5][2]);
          h_bar_r_layer_6->Fill(bar_layer_array[6][2]);
          h_bar_r_layer_7->Fill(bar_layer_array[7][2]);
          h_bar_r_layer_8->Fill(bar_layer_array[8][2]);
          h_bar_r_layer_9->Fill(bar_layer_array[9][2]);
          h_bar_r_layer_10->Fill(bar_layer_array[10][2]);
          h_bar_r_layer_11->Fill(bar_layer_array[11][2]);
          h_bar_r_layer_12->Fill(bar_layer_array[12][2]);
          h_bar_r_layer_13->Fill(bar_layer_array[13][2]);
          h_bar_r_layer_14->Fill(bar_layer_array[14][2]);

	  hit_isMasked->clear();
	  hit_energy->clear();
	  hit_slab->clear();
	  
        } 
	// End of writting histos
	
	// FITS
        // Resolution
        fit_res(h_nhit, mu_nhit[i_energy], sig_nhit[i_energy], res_nhit[i_energy]);
        fit_res(h_sume, mu_sume[i_energy], sig_sume[i_energy], res_sume[i_energy]);
	fit_res(h_weight, mu_weight[i_energy], sig_weight[i_energy], res_weight[i_energy]);
        
	// Moliere fits
	double result_molfit[3] = {0.,0.,0.};
	fit_gauss(h_mol, result_molfit[0], result_molfit[1], result_molfit[2]);
	mu_mol[i_energy] = result_molfit[0];
	mu_error_mol[i_energy] = result_molfit[1];
	sig_mol[i_energy] = 0.5*result_molfit[2];
	
	
	//cout<<"Molière radius: ("<<mu_mol[i_energy]<<" +- "<<sig_mol[i_energy]<<") mm"<<endl;
		
	// Shower profile fits
	// Same method than moliere
	// g_mu, g_sig, l_mpv, l_fwhm
        // nhit
	double nhit_shower_fit_result[N_ECAL_LAYERS][3];
	//cout<<"Fits shower nhit"<<endl;
	
	fit_gauss(h_nhit_layer_0, nhit_shower_fit_result[0][0], nhit_shower_fit_result[0][1], nhit_shower_fit_result[0][2]);
        fit_gauss(h_nhit_layer_1, nhit_shower_fit_result[1][0], nhit_shower_fit_result[1][1], nhit_shower_fit_result[1][2]);
	fit_gauss(h_nhit_layer_2, nhit_shower_fit_result[2][0], nhit_shower_fit_result[2][1], nhit_shower_fit_result[2][2]);
	fit_gauss(h_nhit_layer_3, nhit_shower_fit_result[3][0], nhit_shower_fit_result[3][1], nhit_shower_fit_result[3][2]);
	fit_gauss(h_nhit_layer_4, nhit_shower_fit_result[4][0], nhit_shower_fit_result[4][1], nhit_shower_fit_result[4][2]);
	fit_gauss(h_nhit_layer_5, nhit_shower_fit_result[5][0], nhit_shower_fit_result[5][1], nhit_shower_fit_result[5][2]);
	fit_gauss(h_nhit_layer_6, nhit_shower_fit_result[6][0], nhit_shower_fit_result[6][1], nhit_shower_fit_result[6][2]);
	fit_gauss(h_nhit_layer_7, nhit_shower_fit_result[7][0], nhit_shower_fit_result[7][1], nhit_shower_fit_result[7][2]);
	fit_gauss(h_nhit_layer_8, nhit_shower_fit_result[8][0], nhit_shower_fit_result[8][1], nhit_shower_fit_result[8][2]);
	fit_gauss(h_nhit_layer_9, nhit_shower_fit_result[9][0], nhit_shower_fit_result[9][1], nhit_shower_fit_result[9][2]);
	fit_gauss(h_nhit_layer_10, nhit_shower_fit_result[10][0], nhit_shower_fit_result[10][1], nhit_shower_fit_result[10][2]);
	fit_gauss(h_nhit_layer_11, nhit_shower_fit_result[11][0], nhit_shower_fit_result[11][1], nhit_shower_fit_result[11][2]);
	fit_gauss(h_nhit_layer_12, nhit_shower_fit_result[12][0], nhit_shower_fit_result[12][1], nhit_shower_fit_result[12][2]);
	fit_gauss(h_nhit_layer_13, nhit_shower_fit_result[13][0], nhit_shower_fit_result[13][1], nhit_shower_fit_result[13][2]);
	fit_gauss(h_nhit_layer_14, nhit_shower_fit_result[14][0], nhit_shower_fit_result[14][1], nhit_shower_fit_result[14][2]);
			
	mu_nhit_layer_0[i_energy] = nhit_shower_fit_result[0][0];
	mu_error_nhit_layer_0[i_energy] = nhit_shower_fit_result[0][1];
        sig_nhit_layer_0[i_energy] = 0.5*nhit_shower_fit_result[0][2];
	mu_nhit_layer_1[i_energy] = nhit_shower_fit_result[1][0];
	mu_error_nhit_layer_1[i_energy] = nhit_shower_fit_result[1][1];
        sig_nhit_layer_1[i_energy] = 0.5*nhit_shower_fit_result[1][2];
	mu_nhit_layer_2[i_energy] = nhit_shower_fit_result[2][0];
        mu_error_nhit_layer_2[i_energy] = nhit_shower_fit_result[2][1];
	sig_nhit_layer_2[i_energy] = 0.5*nhit_shower_fit_result[2][2];
	mu_nhit_layer_3[i_energy] = nhit_shower_fit_result[3][0];
        mu_error_nhit_layer_3[i_energy] = nhit_shower_fit_result[3][1];
	sig_nhit_layer_3[i_energy] = 0.5*nhit_shower_fit_result[3][2];
	mu_nhit_layer_4[i_energy] = nhit_shower_fit_result[4][0];
        mu_error_nhit_layer_4[i_energy] = nhit_shower_fit_result[4][1];
	sig_nhit_layer_4[i_energy] = 0.5*nhit_shower_fit_result[4][2];
	mu_nhit_layer_5[i_energy] = nhit_shower_fit_result[5][0];
        mu_error_nhit_layer_5[i_energy] = nhit_shower_fit_result[5][1];
	sig_nhit_layer_5[i_energy] = 0.5*nhit_shower_fit_result[5][2];
	mu_nhit_layer_6[i_energy] = nhit_shower_fit_result[6][0];
        mu_error_nhit_layer_6[i_energy] = nhit_shower_fit_result[6][1];
	sig_nhit_layer_6[i_energy] = 0.5*nhit_shower_fit_result[6][2];
	mu_nhit_layer_7[i_energy] = nhit_shower_fit_result[7][0];
        mu_error_nhit_layer_7[i_energy] = nhit_shower_fit_result[7][1];
	sig_nhit_layer_7[i_energy] = 0.5*nhit_shower_fit_result[7][2];
	mu_nhit_layer_8[i_energy] = nhit_shower_fit_result[8][0];
        mu_error_nhit_layer_8[i_energy] = nhit_shower_fit_result[8][1];
	sig_nhit_layer_8[i_energy] = 0.5*nhit_shower_fit_result[8][2];
	mu_nhit_layer_9[i_energy] = nhit_shower_fit_result[9][0];
        mu_error_nhit_layer_9[i_energy] = nhit_shower_fit_result[9][1];
	sig_nhit_layer_9[i_energy] = 0.5*nhit_shower_fit_result[9][2];
	mu_nhit_layer_10[i_energy] = nhit_shower_fit_result[10][0];
        mu_error_nhit_layer_10[i_energy] = nhit_shower_fit_result[10][1];
	sig_nhit_layer_10[i_energy] = 0.5*nhit_shower_fit_result[10][2];
	mu_nhit_layer_11[i_energy] = nhit_shower_fit_result[11][0];
        mu_error_nhit_layer_11[i_energy] = nhit_shower_fit_result[11][1];
	sig_nhit_layer_11[i_energy] = 0.5*nhit_shower_fit_result[11][2];
	mu_nhit_layer_12[i_energy] = nhit_shower_fit_result[12][0];
        mu_error_nhit_layer_12[i_energy] = nhit_shower_fit_result[12][1];
	sig_nhit_layer_12[i_energy] = 0.5*nhit_shower_fit_result[12][2];
	mu_nhit_layer_13[i_energy] = nhit_shower_fit_result[13][0];
        mu_error_nhit_layer_13[i_energy] = nhit_shower_fit_result[13][1];
	sig_nhit_layer_13[i_energy] = 0.5*nhit_shower_fit_result[13][2];
	mu_nhit_layer_14[i_energy] = nhit_shower_fit_result[14][0];
        mu_error_nhit_layer_14[i_energy] = nhit_shower_fit_result[14][1];
	sig_nhit_layer_14[i_energy] = 0.5*nhit_shower_fit_result[14][2];
    
	// nhit normalized
        double nhit_shower_fit_result_n[N_ECAL_LAYERS][3];
	//cout<<"Fits shower nhit normalized"<<endl;
	
	fit_gauss(h_nhit_layer_n_0, nhit_shower_fit_result_n[0][0], nhit_shower_fit_result_n[0][1], nhit_shower_fit_result_n[0][2]);
        fit_gauss(h_nhit_layer_n_1, nhit_shower_fit_result_n[1][0], nhit_shower_fit_result_n[1][1], nhit_shower_fit_result_n[1][2]);
        fit_gauss(h_nhit_layer_n_2, nhit_shower_fit_result_n[2][0], nhit_shower_fit_result_n[2][1], nhit_shower_fit_result_n[2][2]);
        fit_gauss(h_nhit_layer_n_3, nhit_shower_fit_result_n[3][0], nhit_shower_fit_result_n[3][1], nhit_shower_fit_result_n[3][2]);
        fit_gauss(h_nhit_layer_n_4, nhit_shower_fit_result_n[4][0], nhit_shower_fit_result_n[4][1], nhit_shower_fit_result_n[4][2]);
        fit_gauss(h_nhit_layer_n_5, nhit_shower_fit_result_n[5][0], nhit_shower_fit_result_n[5][1], nhit_shower_fit_result_n[5][2]);
        fit_gauss(h_nhit_layer_n_6, nhit_shower_fit_result_n[6][0], nhit_shower_fit_result_n[6][1], nhit_shower_fit_result_n[6][2]);
        fit_gauss(h_nhit_layer_n_7, nhit_shower_fit_result_n[7][0], nhit_shower_fit_result_n[7][1], nhit_shower_fit_result_n[7][2]);
        fit_gauss(h_nhit_layer_n_8, nhit_shower_fit_result_n[8][0], nhit_shower_fit_result_n[8][1], nhit_shower_fit_result_n[8][2]);
        fit_gauss(h_nhit_layer_n_9, nhit_shower_fit_result_n[9][0], nhit_shower_fit_result_n[9][1], nhit_shower_fit_result_n[9][2]);
        fit_gauss(h_nhit_layer_n_10, nhit_shower_fit_result_n[10][0], nhit_shower_fit_result_n[10][1], nhit_shower_fit_result_n[10][2]);
        fit_gauss(h_nhit_layer_n_11, nhit_shower_fit_result_n[11][0], nhit_shower_fit_result_n[11][1], nhit_shower_fit_result_n[11][2]);
        fit_gauss(h_nhit_layer_n_12, nhit_shower_fit_result_n[12][0], nhit_shower_fit_result_n[12][1], nhit_shower_fit_result_n[12][2]);
        fit_gauss(h_nhit_layer_n_13, nhit_shower_fit_result_n[13][0], nhit_shower_fit_result_n[13][1], nhit_shower_fit_result_n[13][2]);
        fit_gauss(h_nhit_layer_n_14, nhit_shower_fit_result_n[14][0], nhit_shower_fit_result_n[14][1], nhit_shower_fit_result_n[14][2]);
			
        mu_nhit_layer_n_0[i_energy] = nhit_shower_fit_result_n[0][0];
	mu_error_nhit_layer_n_0[i_energy] = nhit_shower_fit_result_n[0][1];
        sig_nhit_layer_n_0[i_energy] = 0.5*nhit_shower_fit_result_n[0][2];
        mu_nhit_layer_n_1[i_energy] = nhit_shower_fit_result_n[1][0];
	mu_error_nhit_layer_n_1[i_energy] = nhit_shower_fit_result_n[1][1];
        sig_nhit_layer_n_1[i_energy] = 0.5*nhit_shower_fit_result_n[1][2];
        mu_nhit_layer_n_2[i_energy] = nhit_shower_fit_result_n[2][0];
	mu_error_nhit_layer_n_2[i_energy] = nhit_shower_fit_result_n[2][1];
        sig_nhit_layer_n_2[i_energy] = 0.5*nhit_shower_fit_result_n[2][2];
        mu_nhit_layer_n_3[i_energy] = nhit_shower_fit_result_n[3][0];
	mu_error_nhit_layer_n_3[i_energy] = nhit_shower_fit_result_n[3][1];
        sig_nhit_layer_n_3[i_energy] = 0.5*nhit_shower_fit_result_n[3][2];
        mu_nhit_layer_n_4[i_energy] = nhit_shower_fit_result_n[4][0];
	mu_error_nhit_layer_n_4[i_energy] = nhit_shower_fit_result_n[4][1];
        sig_nhit_layer_n_4[i_energy] = 0.5*nhit_shower_fit_result_n[4][2];
        mu_nhit_layer_n_5[i_energy] = nhit_shower_fit_result_n[5][0];
	mu_error_nhit_layer_n_5[i_energy] = nhit_shower_fit_result_n[5][1];
        sig_nhit_layer_n_5[i_energy] = 0.5*nhit_shower_fit_result_n[5][2];
        mu_nhit_layer_n_6[i_energy] = nhit_shower_fit_result_n[6][0];
	mu_error_nhit_layer_n_6[i_energy] = nhit_shower_fit_result_n[6][1];
        sig_nhit_layer_n_6[i_energy] = 0.5*nhit_shower_fit_result_n[6][2];
        mu_nhit_layer_n_7[i_energy] = nhit_shower_fit_result_n[7][0];
	mu_error_nhit_layer_n_7[i_energy] = nhit_shower_fit_result_n[7][1];
        sig_nhit_layer_n_7[i_energy] = 0.5*nhit_shower_fit_result_n[7][2];
        mu_nhit_layer_n_8[i_energy] = nhit_shower_fit_result_n[8][0];
	mu_error_nhit_layer_n_8[i_energy] = nhit_shower_fit_result_n[8][1];
        sig_nhit_layer_n_8[i_energy] = 0.5*nhit_shower_fit_result_n[8][2];
        mu_nhit_layer_n_9[i_energy] = nhit_shower_fit_result_n[9][0];
	mu_error_nhit_layer_n_9[i_energy] = nhit_shower_fit_result_n[9][1];
        sig_nhit_layer_n_9[i_energy] = 0.5*nhit_shower_fit_result_n[9][2];
        mu_nhit_layer_n_10[i_energy] = nhit_shower_fit_result_n[10][0];
	mu_error_nhit_layer_n_10[i_energy] = nhit_shower_fit_result_n[10][1];
        sig_nhit_layer_n_10[i_energy] = 0.5*nhit_shower_fit_result_n[10][2];
        mu_nhit_layer_n_11[i_energy] = nhit_shower_fit_result_n[11][0];
	mu_error_nhit_layer_n_11[i_energy] = nhit_shower_fit_result_n[11][1];
        sig_nhit_layer_n_11[i_energy] = 0.5*nhit_shower_fit_result_n[11][2];
        mu_nhit_layer_n_12[i_energy] = nhit_shower_fit_result_n[12][0];
        mu_error_nhit_layer_n_12[i_energy] = nhit_shower_fit_result_n[12][1];
        sig_nhit_layer_n_12[i_energy] = 0.5*nhit_shower_fit_result_n[12][2];
        mu_nhit_layer_n_13[i_energy] = nhit_shower_fit_result_n[13][0];
	mu_error_nhit_layer_n_13[i_energy] = nhit_shower_fit_result_n[13][1];
        sig_nhit_layer_n_13[i_energy] = 0.5*nhit_shower_fit_result_n[13][2];
        mu_nhit_layer_n_14[i_energy] = nhit_shower_fit_result_n[14][0];
	mu_error_nhit_layer_n_14[i_energy] = nhit_shower_fit_result_n[14][1];
        sig_nhit_layer_n_14[i_energy] = 0.5*nhit_shower_fit_result_n[14][2];

	// weight
        double weight_shower_fit_result[N_ECAL_LAYERS][3];
	//cout<<"Fits shower energy weighted"<<endl;
	
	fit_gauss(h_weight_layer_0, weight_shower_fit_result[0][0], weight_shower_fit_result[0][1], weight_shower_fit_result[0][2]);
        fit_gauss(h_weight_layer_1, weight_shower_fit_result[1][0], weight_shower_fit_result[1][1], weight_shower_fit_result[1][2]);
        fit_gauss(h_weight_layer_2, weight_shower_fit_result[2][0], weight_shower_fit_result[2][1], weight_shower_fit_result[2][2]);
        fit_gauss(h_weight_layer_3, weight_shower_fit_result[3][0], weight_shower_fit_result[3][1], weight_shower_fit_result[3][2]);
        fit_gauss(h_weight_layer_4, weight_shower_fit_result[4][0], weight_shower_fit_result[4][1], weight_shower_fit_result[4][2]);
	fit_gauss(h_weight_layer_5, weight_shower_fit_result[5][0], weight_shower_fit_result[5][1], weight_shower_fit_result[5][2]);
        fit_gauss(h_weight_layer_6, weight_shower_fit_result[6][0], weight_shower_fit_result[6][1], weight_shower_fit_result[6][2]);
        fit_gauss(h_weight_layer_7, weight_shower_fit_result[7][0], weight_shower_fit_result[7][1], weight_shower_fit_result[7][2]);
        fit_gauss(h_weight_layer_8, weight_shower_fit_result[8][0], weight_shower_fit_result[8][1], weight_shower_fit_result[8][2]);
        fit_gauss(h_weight_layer_9, weight_shower_fit_result[9][0], weight_shower_fit_result[9][1], weight_shower_fit_result[9][2]);
        fit_gauss(h_weight_layer_10, weight_shower_fit_result[10][0], weight_shower_fit_result[10][1], weight_shower_fit_result[10][2]);
        fit_gauss(h_weight_layer_11, weight_shower_fit_result[11][0], weight_shower_fit_result[11][1], weight_shower_fit_result[11][2]);
        fit_gauss(h_weight_layer_12, weight_shower_fit_result[12][0], weight_shower_fit_result[12][1], weight_shower_fit_result[12][2]);
        fit_gauss(h_weight_layer_13, weight_shower_fit_result[13][0], weight_shower_fit_result[13][1], weight_shower_fit_result[13][2]);
        fit_gauss(h_weight_layer_14, weight_shower_fit_result[14][0], weight_shower_fit_result[14][1], weight_shower_fit_result[14][2]);
	
        mu_weight_layer_0[i_energy] = weight_shower_fit_result[0][0];
        mu_error_weight_layer_0[i_energy] = weight_shower_fit_result[0][1];
        sig_weight_layer_0[i_energy] = 0.5*weight_shower_fit_result[0][2];
        mu_weight_layer_1[i_energy] = weight_shower_fit_result[1][0];
        mu_error_weight_layer_1[i_energy] = weight_shower_fit_result[1][1];
        sig_weight_layer_1[i_energy] = 0.5*weight_shower_fit_result[1][2];
        mu_weight_layer_2[i_energy] = weight_shower_fit_result[2][0];
        mu_error_weight_layer_2[i_energy] = weight_shower_fit_result[2][1];
        sig_weight_layer_2[i_energy] = 0.5*weight_shower_fit_result[2][2];
        mu_weight_layer_3[i_energy] = weight_shower_fit_result[3][0];
        mu_error_weight_layer_3[i_energy] = weight_shower_fit_result[3][1];
        sig_weight_layer_3[i_energy] = 0.5*weight_shower_fit_result[3][2];
        mu_weight_layer_4[i_energy] = weight_shower_fit_result[4][0];
        mu_error_weight_layer_4[i_energy] = weight_shower_fit_result[4][1];
	sig_weight_layer_4[i_energy] = 0.5*weight_shower_fit_result[4][2];
        mu_weight_layer_5[i_energy] = weight_shower_fit_result[5][0];
        mu_error_weight_layer_5[i_energy] = weight_shower_fit_result[5][1];
        sig_weight_layer_5[i_energy] = 0.5*weight_shower_fit_result[5][2];
        mu_weight_layer_6[i_energy] = weight_shower_fit_result[6][0];
        mu_error_weight_layer_6[i_energy] = weight_shower_fit_result[6][1];
        sig_weight_layer_6[i_energy] = 0.5*weight_shower_fit_result[6][2];
        mu_weight_layer_7[i_energy] = weight_shower_fit_result[7][0];
        mu_error_weight_layer_7[i_energy] = weight_shower_fit_result[7][1];
        sig_weight_layer_7[i_energy] = 0.5*weight_shower_fit_result[7][2];
        mu_weight_layer_8[i_energy] = weight_shower_fit_result[8][0];
        mu_error_weight_layer_8[i_energy] = weight_shower_fit_result[8][1];
        sig_weight_layer_8[i_energy] = 0.5*weight_shower_fit_result[8][2];
        mu_weight_layer_9[i_energy] = weight_shower_fit_result[9][0];
        mu_error_weight_layer_9[i_energy] = weight_shower_fit_result[9][1];
        sig_weight_layer_9[i_energy] = 0.5*weight_shower_fit_result[9][2];
        mu_weight_layer_10[i_energy] = weight_shower_fit_result[10][0];
        mu_error_weight_layer_10[i_energy] = weight_shower_fit_result[10][1];
        sig_weight_layer_10[i_energy] = 0.5*weight_shower_fit_result[10][2];
        mu_weight_layer_11[i_energy] = weight_shower_fit_result[11][0];
        mu_error_weight_layer_11[i_energy] = weight_shower_fit_result[11][1];
        sig_weight_layer_11[i_energy] = 0.5*weight_shower_fit_result[11][2];
        mu_weight_layer_12[i_energy] = weight_shower_fit_result[12][0];
        mu_error_weight_layer_12[i_energy] = weight_shower_fit_result[12][1];
        sig_weight_layer_12[i_energy] = 0.5*weight_shower_fit_result[12][2];
        mu_weight_layer_13[i_energy] = weight_shower_fit_result[13][0];
        mu_error_weight_layer_13[i_energy] = weight_shower_fit_result[13][1];
        sig_weight_layer_13[i_energy] = 0.5*weight_shower_fit_result[13][2];
        mu_weight_layer_14[i_energy] = weight_shower_fit_result[14][0];
        mu_error_weight_layer_14[i_energy] = weight_shower_fit_result[14][1];
        sig_weight_layer_14[i_energy] = 0.5*weight_shower_fit_result[14][2];

	// weight normalized
        double weight_shower_fit_result_n[N_ECAL_LAYERS][3];
        //cout<<"Fits shower energy weighted normalized"<<endl;
	
	fit_gauss(h_weight_layer_n_0, weight_shower_fit_result_n[0][0], weight_shower_fit_result_n[0][1], weight_shower_fit_result_n[0][2]);
        fit_gauss(h_weight_layer_n_1, weight_shower_fit_result_n[1][0], weight_shower_fit_result_n[1][1], weight_shower_fit_result_n[1][2]);
        fit_gauss(h_weight_layer_n_2, weight_shower_fit_result_n[2][0], weight_shower_fit_result_n[2][1], weight_shower_fit_result_n[2][2]);
        fit_gauss(h_weight_layer_n_3, weight_shower_fit_result_n[3][0], weight_shower_fit_result_n[3][1], weight_shower_fit_result_n[3][2]);
        fit_gauss(h_weight_layer_n_4, weight_shower_fit_result_n[4][0], weight_shower_fit_result_n[4][1], weight_shower_fit_result_n[4][2]);
        fit_gauss(h_weight_layer_n_5, weight_shower_fit_result_n[5][0], weight_shower_fit_result_n[5][1], weight_shower_fit_result_n[5][2]);
        fit_gauss(h_weight_layer_n_6, weight_shower_fit_result_n[6][0], weight_shower_fit_result_n[6][1], weight_shower_fit_result_n[6][2]);
        fit_gauss(h_weight_layer_n_7, weight_shower_fit_result_n[7][0], weight_shower_fit_result_n[7][1], weight_shower_fit_result_n[7][2]);
        fit_gauss(h_weight_layer_n_8, weight_shower_fit_result_n[8][0], weight_shower_fit_result_n[8][1], weight_shower_fit_result_n[8][2]);
        fit_gauss(h_weight_layer_n_9, weight_shower_fit_result_n[9][0], weight_shower_fit_result_n[9][1], weight_shower_fit_result_n[9][2]);
        fit_gauss(h_weight_layer_n_10, weight_shower_fit_result_n[10][0], weight_shower_fit_result_n[10][1], weight_shower_fit_result_n[10][2]);
        fit_gauss(h_weight_layer_n_11, weight_shower_fit_result_n[11][0], weight_shower_fit_result_n[11][1], weight_shower_fit_result_n[11][2]);
        fit_gauss(h_weight_layer_n_12, weight_shower_fit_result_n[12][0], weight_shower_fit_result_n[12][1], weight_shower_fit_result_n[12][2]);
        fit_gauss(h_weight_layer_n_13, weight_shower_fit_result_n[13][0], weight_shower_fit_result_n[13][1], weight_shower_fit_result_n[13][2]);
        fit_gauss(h_weight_layer_n_14, weight_shower_fit_result_n[14][0], weight_shower_fit_result_n[14][1], weight_shower_fit_result_n[14][2]);
	
        mu_weight_layer_n_0[i_energy] = weight_shower_fit_result_n[0][0];
        mu_error_weight_layer_n_0[i_energy] = weight_shower_fit_result_n[0][1];
        sig_weight_layer_n_0[i_energy] = 0.5*weight_shower_fit_result_n[0][2];
        mu_weight_layer_n_1[i_energy] = weight_shower_fit_result_n[1][0];
        mu_error_weight_layer_n_1[i_energy] = weight_shower_fit_result_n[1][1];
        sig_weight_layer_n_1[i_energy] = 0.5*weight_shower_fit_result_n[1][2];
        mu_weight_layer_n_2[i_energy] = weight_shower_fit_result_n[2][0];
        mu_error_weight_layer_n_2[i_energy] = weight_shower_fit_result_n[2][1];
        sig_weight_layer_n_2[i_energy] = 0.5*weight_shower_fit_result_n[2][2];
        mu_weight_layer_n_3[i_energy] = weight_shower_fit_result_n[3][0];
        mu_error_weight_layer_n_3[i_energy] = weight_shower_fit_result_n[3][1];
        sig_weight_layer_n_3[i_energy] = 0.5*weight_shower_fit_result_n[3][2];
        mu_weight_layer_n_4[i_energy] = weight_shower_fit_result_n[4][0];
        mu_error_weight_layer_n_4[i_energy] = weight_shower_fit_result_n[4][1];
        sig_weight_layer_n_4[i_energy] = 0.5*weight_shower_fit_result_n[4][2];
        mu_weight_layer_n_5[i_energy] = weight_shower_fit_result_n[5][0];
        mu_error_weight_layer_n_5[i_energy] = weight_shower_fit_result_n[5][1];
        sig_weight_layer_n_5[i_energy] = 0.5*weight_shower_fit_result_n[5][2];
        mu_weight_layer_n_6[i_energy] = weight_shower_fit_result_n[6][0];
        mu_error_weight_layer_n_6[i_energy] = weight_shower_fit_result_n[6][1];
        sig_weight_layer_n_6[i_energy] = 0.5*weight_shower_fit_result_n[6][2];
        mu_weight_layer_n_7[i_energy] = weight_shower_fit_result_n[7][0];
        mu_error_weight_layer_n_7[i_energy] = weight_shower_fit_result_n[7][1];
        sig_weight_layer_n_7[i_energy] = 0.5*weight_shower_fit_result_n[7][2];
        mu_weight_layer_n_8[i_energy] = weight_shower_fit_result_n[8][0];
        mu_error_weight_layer_n_8[i_energy] = weight_shower_fit_result_n[8][1];
        sig_weight_layer_n_8[i_energy] = 0.5*weight_shower_fit_result_n[8][2];
        mu_weight_layer_n_9[i_energy] = weight_shower_fit_result_n[9][0];
        mu_error_weight_layer_n_9[i_energy] = weight_shower_fit_result_n[9][1];
        sig_weight_layer_n_9[i_energy] = 0.5*weight_shower_fit_result_n[9][2];
        mu_weight_layer_n_10[i_energy] = weight_shower_fit_result_n[10][0];
        mu_error_weight_layer_n_10[i_energy] = weight_shower_fit_result_n[10][1];
        sig_weight_layer_n_10[i_energy] = 0.5*weight_shower_fit_result_n[10][2];
        mu_weight_layer_n_11[i_energy] = weight_shower_fit_result_n[11][0];
        mu_error_weight_layer_n_11[i_energy] = weight_shower_fit_result_n[11][1];
        sig_weight_layer_n_11[i_energy] = 0.5*weight_shower_fit_result_n[11][2];
        mu_weight_layer_n_12[i_energy] = weight_shower_fit_result_n[12][0];
        mu_error_weight_layer_n_12[i_energy] = weight_shower_fit_result_n[12][1];
        sig_weight_layer_n_12[i_energy] = 0.5*weight_shower_fit_result_n[12][2];
        mu_weight_layer_n_13[i_energy] = weight_shower_fit_result_n[13][0];
        mu_error_weight_layer_n_13[i_energy] = weight_shower_fit_result_n[13][1];
        sig_weight_layer_n_13[i_energy] = 0.5*weight_shower_fit_result_n[13][2];
        mu_weight_layer_n_14[i_energy] = weight_shower_fit_result_n[14][0];
        mu_error_weight_layer_n_14[i_energy] = weight_shower_fit_result_n[14][1];
        sig_weight_layer_n_14[i_energy] = 0.5*weight_shower_fit_result_n[14][2];
	
	// barycenter x per layer
        double bar_x_fit_result[N_ECAL_LAYERS][2];
        
	fit_bar(h_bar_x_layer_0, bar_x_fit_result[0][0], bar_x_fit_result[0][1]);
	fit_bar(h_bar_x_layer_1, bar_x_fit_result[1][0], bar_x_fit_result[1][1]);
	fit_bar(h_bar_x_layer_2, bar_x_fit_result[2][0], bar_x_fit_result[2][1]);
	fit_bar(h_bar_x_layer_3, bar_x_fit_result[3][0], bar_x_fit_result[3][1]);
	fit_bar(h_bar_x_layer_4, bar_x_fit_result[4][0], bar_x_fit_result[4][1]);
	fit_bar(h_bar_x_layer_5, bar_x_fit_result[5][0], bar_x_fit_result[5][1]);
	fit_bar(h_bar_x_layer_6, bar_x_fit_result[6][0], bar_x_fit_result[6][1]);
	fit_bar(h_bar_x_layer_7, bar_x_fit_result[7][0], bar_x_fit_result[7][1]);
	fit_bar(h_bar_x_layer_8, bar_x_fit_result[8][0], bar_x_fit_result[8][1]);
	fit_bar(h_bar_x_layer_9, bar_x_fit_result[9][0], bar_x_fit_result[9][1]);
	fit_bar(h_bar_x_layer_10, bar_x_fit_result[10][0], bar_x_fit_result[10][1]);
	fit_bar(h_bar_x_layer_11, bar_x_fit_result[11][0], bar_x_fit_result[11][1]);
	fit_bar(h_bar_x_layer_12, bar_x_fit_result[12][0], bar_x_fit_result[12][1]);
	fit_bar(h_bar_x_layer_13, bar_x_fit_result[13][0], bar_x_fit_result[13][1]);
	fit_bar(h_bar_x_layer_14, bar_x_fit_result[14][0], bar_x_fit_result[14][1]);
	

        mu_barycenter_x_layer_0[i_energy] = bar_x_fit_result[0][0];
        sig_barycenter_x_layer_0[i_energy] = bar_x_fit_result[0][1];
	mu_barycenter_x_layer_1[i_energy] = bar_x_fit_result[1][0];
        sig_barycenter_x_layer_1[i_energy] = bar_x_fit_result[1][1];
	mu_barycenter_x_layer_2[i_energy] = bar_x_fit_result[2][0];
        sig_barycenter_x_layer_2[i_energy] = bar_x_fit_result[2][1];
	mu_barycenter_x_layer_3[i_energy] = bar_x_fit_result[3][0];
        sig_barycenter_x_layer_3[i_energy] = bar_x_fit_result[3][1];
	mu_barycenter_x_layer_4[i_energy] = bar_x_fit_result[4][0];
        sig_barycenter_x_layer_4[i_energy] = bar_x_fit_result[4][1];
	mu_barycenter_x_layer_5[i_energy] = bar_x_fit_result[5][0];
        sig_barycenter_x_layer_5[i_energy] = bar_x_fit_result[5][1];
	mu_barycenter_x_layer_6[i_energy] = bar_x_fit_result[6][0];
        sig_barycenter_x_layer_6[i_energy] = bar_x_fit_result[6][1];
	mu_barycenter_x_layer_7[i_energy] = bar_x_fit_result[7][0];
        sig_barycenter_x_layer_7[i_energy] = bar_x_fit_result[7][1];
	mu_barycenter_x_layer_8[i_energy] = bar_x_fit_result[8][0];
        sig_barycenter_x_layer_8[i_energy] = bar_x_fit_result[8][1];
	mu_barycenter_x_layer_9[i_energy] = bar_x_fit_result[9][0];
        sig_barycenter_x_layer_9[i_energy] = bar_x_fit_result[9][1];
	mu_barycenter_x_layer_10[i_energy] = bar_x_fit_result[10][0];
        sig_barycenter_x_layer_10[i_energy] = bar_x_fit_result[10][1];
	mu_barycenter_x_layer_11[i_energy] = bar_x_fit_result[11][0];
        sig_barycenter_x_layer_11[i_energy] = bar_x_fit_result[11][1];
	mu_barycenter_x_layer_12[i_energy] = bar_x_fit_result[12][0];
        sig_barycenter_x_layer_12[i_energy] = bar_x_fit_result[12][1];
	mu_barycenter_x_layer_13[i_energy] = bar_x_fit_result[13][0];
        sig_barycenter_x_layer_13[i_energy] = bar_x_fit_result[13][1];
	mu_barycenter_x_layer_14[i_energy] = bar_x_fit_result[14][0];
        sig_barycenter_x_layer_14[i_energy] = bar_x_fit_result[14][1];

	// barycenter y per layer 
        double bar_y_fit_result[N_ECAL_LAYERS][2];
        
        fit_bar(h_bar_y_layer_0, bar_y_fit_result[0][0], bar_y_fit_result[0][1]);
        fit_bar(h_bar_y_layer_1, bar_y_fit_result[1][0], bar_y_fit_result[1][1]);
        fit_bar(h_bar_y_layer_2, bar_y_fit_result[2][0], bar_y_fit_result[2][1]);
        fit_bar(h_bar_y_layer_3, bar_y_fit_result[3][0], bar_y_fit_result[3][1]);
        fit_bar(h_bar_y_layer_4, bar_y_fit_result[4][0], bar_y_fit_result[4][1]);
        fit_bar(h_bar_y_layer_5, bar_y_fit_result[5][0], bar_y_fit_result[5][1]);
        fit_bar(h_bar_y_layer_6, bar_y_fit_result[6][0], bar_y_fit_result[6][1]);
        fit_bar(h_bar_y_layer_7, bar_y_fit_result[7][0], bar_y_fit_result[7][1]);
        fit_bar(h_bar_y_layer_8, bar_y_fit_result[8][0], bar_y_fit_result[8][1]);
        fit_bar(h_bar_y_layer_9, bar_y_fit_result[9][0], bar_y_fit_result[9][1]);
        fit_bar(h_bar_y_layer_10, bar_y_fit_result[10][0], bar_y_fit_result[10][1]);
        fit_bar(h_bar_y_layer_11, bar_y_fit_result[11][0], bar_y_fit_result[11][1]);
        fit_bar(h_bar_y_layer_12, bar_y_fit_result[12][0], bar_y_fit_result[12][1]);
        fit_bar(h_bar_y_layer_13, bar_y_fit_result[13][0], bar_y_fit_result[13][1]);
        fit_bar(h_bar_y_layer_14, bar_y_fit_result[14][0], bar_y_fit_result[14][1]);
	

        mu_barycenter_y_layer_0[i_energy] = bar_y_fit_result[0][0];
        sig_barycenter_y_layer_0[i_energy] = bar_y_fit_result[0][1];
        mu_barycenter_y_layer_1[i_energy] = bar_y_fit_result[1][0];
        sig_barycenter_y_layer_1[i_energy] = bar_y_fit_result[1][1];
        mu_barycenter_y_layer_2[i_energy] = bar_y_fit_result[2][0];
        sig_barycenter_y_layer_2[i_energy] = bar_y_fit_result[2][1];
        mu_barycenter_y_layer_3[i_energy] = bar_y_fit_result[3][0];
        sig_barycenter_y_layer_3[i_energy] = bar_y_fit_result[3][1];
        mu_barycenter_y_layer_4[i_energy] = bar_y_fit_result[4][0];
        sig_barycenter_y_layer_4[i_energy] = bar_y_fit_result[4][1];
        mu_barycenter_y_layer_5[i_energy] = bar_y_fit_result[5][0];
        sig_barycenter_y_layer_5[i_energy] = bar_y_fit_result[5][1];
        mu_barycenter_y_layer_6[i_energy] = bar_y_fit_result[6][0];
        sig_barycenter_y_layer_6[i_energy] = bar_y_fit_result[6][1];
        mu_barycenter_y_layer_7[i_energy] = bar_y_fit_result[7][0];
        sig_barycenter_y_layer_7[i_energy] = bar_y_fit_result[7][1];
        mu_barycenter_y_layer_8[i_energy] = bar_y_fit_result[8][0];
        sig_barycenter_y_layer_8[i_energy] = bar_y_fit_result[8][1];
        mu_barycenter_y_layer_9[i_energy] = bar_y_fit_result[9][0];
        sig_barycenter_y_layer_9[i_energy] = bar_y_fit_result[9][1];
        mu_barycenter_y_layer_10[i_energy] = bar_y_fit_result[10][0];
        sig_barycenter_y_layer_10[i_energy] = bar_y_fit_result[10][1];
        mu_barycenter_y_layer_11[i_energy] = bar_y_fit_result[11][0];
        sig_barycenter_y_layer_11[i_energy] = bar_y_fit_result[11][1];
        mu_barycenter_y_layer_12[i_energy] = bar_y_fit_result[12][0];
        sig_barycenter_y_layer_12[i_energy] = bar_y_fit_result[12][1];
        mu_barycenter_y_layer_13[i_energy] = bar_y_fit_result[13][0];
        sig_barycenter_y_layer_13[i_energy] = bar_y_fit_result[13][1];
        mu_barycenter_y_layer_14[i_energy] = bar_y_fit_result[14][0];
        sig_barycenter_y_layer_14[i_energy] = bar_y_fit_result[14][1];

	// barycenter r per layer
        double bar_r_fit_result[N_ECAL_LAYERS][2];

        fit_bar(h_bar_r_layer_0, bar_r_fit_result[0][0], bar_r_fit_result[0][1]);
        fit_bar(h_bar_r_layer_1, bar_r_fit_result[1][0], bar_r_fit_result[1][1]);
        fit_bar(h_bar_r_layer_2, bar_r_fit_result[2][0], bar_r_fit_result[2][1]);
        fit_bar(h_bar_r_layer_3, bar_r_fit_result[3][0], bar_r_fit_result[3][1]);
        fit_bar(h_bar_r_layer_4, bar_r_fit_result[4][0], bar_r_fit_result[4][1]);
        fit_bar(h_bar_r_layer_5, bar_r_fit_result[5][0], bar_r_fit_result[5][1]);
        fit_bar(h_bar_r_layer_6, bar_r_fit_result[6][0], bar_r_fit_result[6][1]);
        fit_bar(h_bar_r_layer_7, bar_r_fit_result[7][0], bar_r_fit_result[7][1]);
        fit_bar(h_bar_r_layer_8, bar_r_fit_result[8][0], bar_r_fit_result[8][1]);
        fit_bar(h_bar_r_layer_9, bar_r_fit_result[9][0], bar_r_fit_result[9][1]);
        fit_bar(h_bar_r_layer_10, bar_r_fit_result[10][0], bar_r_fit_result[10][1]);
        fit_bar(h_bar_r_layer_11, bar_r_fit_result[11][0], bar_r_fit_result[11][1]);
        fit_bar(h_bar_r_layer_12, bar_r_fit_result[12][0], bar_r_fit_result[12][1]);
        fit_bar(h_bar_r_layer_13, bar_r_fit_result[13][0], bar_r_fit_result[13][1]);
        fit_bar(h_bar_r_layer_14, bar_r_fit_result[14][0], bar_r_fit_result[14][1]);


        mu_barycenter_r_layer_0[i_energy] = bar_r_fit_result[0][0];
        sig_barycenter_r_layer_0[i_energy] = bar_r_fit_result[0][1];
        mu_barycenter_r_layer_1[i_energy] = bar_r_fit_result[1][0];
        sig_barycenter_r_layer_1[i_energy] = bar_r_fit_result[1][1];
        mu_barycenter_r_layer_2[i_energy] = bar_r_fit_result[2][0];
        sig_barycenter_r_layer_2[i_energy] = bar_r_fit_result[2][1];
        mu_barycenter_r_layer_3[i_energy] = bar_r_fit_result[3][0];
        sig_barycenter_r_layer_3[i_energy] = bar_r_fit_result[3][1];
        mu_barycenter_r_layer_4[i_energy] = bar_r_fit_result[4][0];
        sig_barycenter_r_layer_4[i_energy] = bar_r_fit_result[4][1];
        mu_barycenter_r_layer_5[i_energy] = bar_r_fit_result[5][0];
        sig_barycenter_r_layer_5[i_energy] = bar_r_fit_result[5][1];
        mu_barycenter_r_layer_6[i_energy] = bar_r_fit_result[6][0];
        sig_barycenter_r_layer_6[i_energy] = bar_r_fit_result[6][1];
        mu_barycenter_r_layer_7[i_energy] = bar_r_fit_result[7][0];
        sig_barycenter_r_layer_7[i_energy] = bar_r_fit_result[7][1];
        mu_barycenter_r_layer_8[i_energy] = bar_r_fit_result[8][0];
        sig_barycenter_r_layer_8[i_energy] = bar_r_fit_result[8][1];
        mu_barycenter_r_layer_9[i_energy] = bar_r_fit_result[9][0];
        sig_barycenter_r_layer_9[i_energy] = bar_r_fit_result[9][1];
        mu_barycenter_r_layer_10[i_energy] = bar_r_fit_result[10][0];
        sig_barycenter_r_layer_10[i_energy] = bar_r_fit_result[10][1];
        mu_barycenter_r_layer_11[i_energy] = bar_r_fit_result[11][0];
        sig_barycenter_r_layer_11[i_energy] = bar_r_fit_result[11][1];
        mu_barycenter_r_layer_12[i_energy] = bar_r_fit_result[12][0];
        sig_barycenter_r_layer_12[i_energy] = bar_r_fit_result[12][1];
        mu_barycenter_r_layer_13[i_energy] = bar_r_fit_result[13][0];
        sig_barycenter_r_layer_13[i_energy] = bar_r_fit_result[13][1];
        mu_barycenter_r_layer_14[i_energy] = bar_r_fit_result[14][0];
        sig_barycenter_r_layer_14[i_energy] = bar_r_fit_result[14][1];

	//cout<<"barycenter x laxer: "<<bar_x_fit_result[0][0]<<", "<<bar_x_fit_result[1][0]<<", "<<bar_x_fit_result[2][0]<<", "<<bar_x_fit_result[3][0]<<", "<<bar_x_fit_result[4][0]<<", "<<bar_x_fit_result[5][0]<<", "<<bar_x_fit_result[6][0]<<", "<<bar_x_fit_result[7][0]<<", "<<bar_x_fit_result[8][0]<<", "<<bar_x_fit_result[9][0]<<", "<<bar_x_fit_result[10][0]<<", "<<bar_x_fit_result[11][0]<<", "<<bar_x_fit_result[12][0]<<", "<<bar_x_fit_result[13][0]<<", "<<bar_x_fit_result[14][0]<<", "<<endl;

	//cout<<"barycenter y layer: "<<bar_y_fit_result[0][0]<<", "<<bar_y_fit_result[1][0]<<", "<<bar_y_fit_result[2][0]<<", "<<bar_y_fit_result[3][0]<<", "<<bar_y_fit_result[4][0]<<", "<<bar_y_fit_result[5][0]<<", "<<bar_y_fit_result[6][0]<<", "<<bar_y_fit_result[7][0]<<", "<<bar_y_fit_result[8][0]<<", "<<bar_y_fit_result[9][0]<<", "<<bar_y_fit_result[10][0]<<", "<<bar_y_fit_result[11][0]<<", "<<bar_y_fit_result[12][0]<<", "<<bar_y_fit_result[13][0]<<", "<<bar_y_fit_result[14][0]<<", "<<endl;

	// The rest of variables
	// If use the maximum bin remember -1 because the binning start at 1!
	// or -2 for the layers since bin 1 is layer -1
	
	shower_nhit_max_layer[i_energy] = h_shower_nhit_max_layer->GetMaximumBin()-2;
	shower_nhit_start_layer[i_energy] = h_shower_nhit_start_layer->GetMaximumBin()-2;
	shower_nhit_end_layer[i_energy] = h_shower_nhit_end_layer->GetMaximumBin()-2;
	shower_nhit_start_10_layer[i_energy] = h_shower_nhit_start_10_layer->GetMaximumBin()-2;
        shower_nhit_end_10_layer[i_energy] = h_shower_nhit_end_10_layer->GetMaximumBin()-2;
	shower_nhit_average[i_energy] = h_shower_nhit_average->GetMaximumBin()-1; 
	shower_nhit_max[i_energy] = h_shower_nhit_max->GetMaximumBin()-1;
	
	//cout<<"Nhit shower layer start: "<<shower_nhit_start_layer[i_energy]<<", max: "<<shower_nhit_max_layer[i_energy]<<", end: "<<shower_nhit_end_layer[i_energy]<<endl;
	//cout<<"Nhit shower layer start (0.1*Max): "<<shower_nhit_start_10_layer[i_energy]<<", end (0.1*Max): "<<shower_nhit_end_10_layer[i_energy]<<endl;
	//cout<<"Nhit shower layer average: "<<shower_nhit_average[i_energy]<<", max. value: "<<shower_nhit_max[i_energy]<<endl;
	
	shower_weight_max_layer[i_energy] = h_shower_weight_max_layer->GetMaximumBin()-2;
        shower_weight_start_layer[i_energy] = h_shower_weight_start_layer->GetMaximumBin()-2;
        shower_weight_end_layer[i_energy] = h_shower_weight_end_layer->GetMaximumBin()-2;
	shower_weight_start_10_layer[i_energy] = h_shower_weight_start_10_layer->GetMaximumBin()-2;
        shower_weight_end_10_layer[i_energy] = h_shower_weight_end_10_layer->GetMaximumBin()-2;
        shower_weight_average[i_energy] = h_shower_weight_average->GetMaximumBin()-1;
        shower_weight_max[i_energy] = h_shower_weight_max->GetMaximumBin()-1;

        //cout<<"Weight shower layer start: "<<shower_weight_start_layer[i_energy]<<", max: "<<shower_weight_max_layer[i_energy]<<", end: "<<shower_weight_end_layer[i_energy]<<endl;
	//cout<<"Weight shower layer start (0.1*Max): "<<shower_weight_start_10_layer[i_energy]<<", end (0.1*Max): "<<shower_weight_end_10_layer[i_energy]<<endl;
        //cout<<"Weight shower layer average: "<<shower_weight_average[i_energy]<<", max. value: "<<shower_weight_max[i_energy]<<endl;
	
	double bar_fits[4][2];

        fit_bar(h_bar_x, bar_fits[0][0], bar_fits[0][1]);
	fit_bar(h_bar_y, bar_fits[1][0], bar_fits[1][1]);
	fit_bar(h_bar_z, bar_fits[2][0], bar_fits[2][1]);
	fit_bar(h_bar_r, bar_fits[3][0], bar_fits[3][1]);

	mu_barycenter_x[i_energy] = bar_fits[0][0];
	mu_barycenter_y[i_energy] = bar_fits[1][0];
	mu_barycenter_z[i_energy] = bar_fits[2][0];
	mu_barycenter_r[i_energy] = bar_fits[3][0];
	sig_barycenter_x[i_energy] = bar_fits[0][1];
        sig_barycenter_y[i_energy] = bar_fits[1][1];
        sig_barycenter_z[i_energy] = bar_fits[2][1];
	sig_barycenter_r[i_energy] = bar_fits[2][1];
	
	//Write objects
        f.WriteTObject(h_nhit);    
        f.WriteTObject(h_sume);    
        f.WriteTObject(h_weight);  
        f.WriteTObject(h_mol);     
	
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
	
	f.WriteTObject(h_nhit_layer_n_0);
        f.WriteTObject(h_nhit_layer_n_1);
        f.WriteTObject(h_nhit_layer_n_2);
        f.WriteTObject(h_nhit_layer_n_3);
        f.WriteTObject(h_nhit_layer_n_4);
        f.WriteTObject(h_nhit_layer_n_5);
        f.WriteTObject(h_nhit_layer_n_6);
        f.WriteTObject(h_nhit_layer_n_7);
        f.WriteTObject(h_nhit_layer_n_8);
        f.WriteTObject(h_nhit_layer_n_9);
        f.WriteTObject(h_nhit_layer_n_10);
        f.WriteTObject(h_nhit_layer_n_11);
        f.WriteTObject(h_nhit_layer_n_12);
        f.WriteTObject(h_nhit_layer_n_13);
        f.WriteTObject(h_nhit_layer_n_14);

	f.WriteTObject(h_shower_nhit_max_layer);
        f.WriteTObject(h_shower_nhit_start_layer);
        f.WriteTObject(h_shower_nhit_end_layer);
        f.WriteTObject(h_shower_nhit_start_10_layer);
        f.WriteTObject(h_shower_nhit_end_10_layer);
	f.WriteTObject(h_shower_nhit_average);
        f.WriteTObject(h_shower_nhit_max);
	
	f.WriteTObject(h_weight_layer_0);
        f.WriteTObject(h_weight_layer_1);
        f.WriteTObject(h_weight_layer_2);
        f.WriteTObject(h_weight_layer_3);
        f.WriteTObject(h_weight_layer_4);
        f.WriteTObject(h_weight_layer_5);
        f.WriteTObject(h_weight_layer_6);
        f.WriteTObject(h_weight_layer_7);
        f.WriteTObject(h_weight_layer_8);
        f.WriteTObject(h_weight_layer_9);
        f.WriteTObject(h_weight_layer_10);
        f.WriteTObject(h_weight_layer_11);
        f.WriteTObject(h_weight_layer_12);
        f.WriteTObject(h_weight_layer_13);
        f.WriteTObject(h_weight_layer_14);

        f.WriteTObject(h_weight_layer_n_0);
        f.WriteTObject(h_weight_layer_n_1);
        f.WriteTObject(h_weight_layer_n_2);
        f.WriteTObject(h_weight_layer_n_3);
        f.WriteTObject(h_weight_layer_n_4);
        f.WriteTObject(h_weight_layer_n_5);
        f.WriteTObject(h_weight_layer_n_6);
        f.WriteTObject(h_weight_layer_n_7);
        f.WriteTObject(h_weight_layer_n_8);
        f.WriteTObject(h_weight_layer_n_9);
        f.WriteTObject(h_weight_layer_n_10);
        f.WriteTObject(h_weight_layer_n_11);
        f.WriteTObject(h_weight_layer_n_12);
        f.WriteTObject(h_weight_layer_n_13);
        f.WriteTObject(h_weight_layer_n_14);

        f.WriteTObject(h_shower_weight_max_layer);
        f.WriteTObject(h_shower_weight_start_layer);
        f.WriteTObject(h_shower_weight_end_layer);
	f.WriteTObject(h_shower_weight_start_10_layer);
        f.WriteTObject(h_shower_weight_end_10_layer);
        f.WriteTObject(h_shower_weight_average);
        f.WriteTObject(h_shower_weight_max);

	f.WriteTObject(h_bar_x_layer_0);
        f.WriteTObject(h_bar_x_layer_1);
        f.WriteTObject(h_bar_x_layer_2);
        f.WriteTObject(h_bar_x_layer_3);
        f.WriteTObject(h_bar_x_layer_4);
        f.WriteTObject(h_bar_x_layer_5);
        f.WriteTObject(h_bar_x_layer_6);
        f.WriteTObject(h_bar_x_layer_7);
        f.WriteTObject(h_bar_x_layer_8);
        f.WriteTObject(h_bar_x_layer_9);
        f.WriteTObject(h_bar_x_layer_10);
        f.WriteTObject(h_bar_x_layer_11);
        f.WriteTObject(h_bar_x_layer_12);
        f.WriteTObject(h_bar_x_layer_13);
        f.WriteTObject(h_bar_x_layer_14);

	f.WriteTObject(h_bar_y_layer_0);
        f.WriteTObject(h_bar_y_layer_1);
        f.WriteTObject(h_bar_y_layer_2);
        f.WriteTObject(h_bar_y_layer_3);
        f.WriteTObject(h_bar_y_layer_4);
        f.WriteTObject(h_bar_y_layer_5);
        f.WriteTObject(h_bar_y_layer_6);
        f.WriteTObject(h_bar_y_layer_7);
        f.WriteTObject(h_bar_y_layer_8);
        f.WriteTObject(h_bar_y_layer_9);
        f.WriteTObject(h_bar_y_layer_10);
        f.WriteTObject(h_bar_y_layer_11);
        f.WriteTObject(h_bar_y_layer_12);
        f.WriteTObject(h_bar_y_layer_13);
        f.WriteTObject(h_bar_y_layer_14);

	f.WriteTObject(h_bar_r_layer_0);
        f.WriteTObject(h_bar_r_layer_1);
        f.WriteTObject(h_bar_r_layer_2);
        f.WriteTObject(h_bar_r_layer_3);
        f.WriteTObject(h_bar_r_layer_4);
        f.WriteTObject(h_bar_r_layer_5);
        f.WriteTObject(h_bar_r_layer_6);
        f.WriteTObject(h_bar_r_layer_7);
        f.WriteTObject(h_bar_r_layer_8);
        f.WriteTObject(h_bar_r_layer_9);
        f.WriteTObject(h_bar_r_layer_10);
        f.WriteTObject(h_bar_r_layer_11);
        f.WriteTObject(h_bar_r_layer_12);
        f.WriteTObject(h_bar_r_layer_13);
        f.WriteTObject(h_bar_r_layer_14);

	f.WriteTObject(h_bar_x);
	f.WriteTObject(h_bar_y);
	f.WriteTObject(h_bar_z);
	f.WriteTObject(h_bar_r);
    }
    
    f.WriteObject(&mu_nhit, "mu_nhit");                     f.WriteObject(&sig_nhit, "sig_nhit");                   f.WriteObject(&res_nhit, "res_nhit");
    f.WriteObject(&mu_sume, "mu_sume");                     f.WriteObject(&sig_sume, "sig_sume");                   f.WriteObject(&res_sume, "res_sume");
    f.WriteObject(&mu_weight, "mu_weight");                 f.WriteObject(&sig_weight, "sig_weight");               f.WriteObject(&res_weight, "res_weight");
    f.WriteObject(&mu_mol, "mu_mol");  f.WriteObject(&sig_mol, "sig_mol");
    
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

    f.WriteTObject(&mu_nhit_layer_n_0, "mu_nhit_layer_n_0");
    f.WriteTObject(&mu_nhit_layer_n_1, "mu_nhit_layer_n_1");
    f.WriteTObject(&mu_nhit_layer_n_2, "mu_nhit_layer_n_2");
    f.WriteTObject(&mu_nhit_layer_n_3, "mu_nhit_layer_n_3");
    f.WriteTObject(&mu_nhit_layer_n_4, "mu_nhit_layer_n_4");
    f.WriteTObject(&mu_nhit_layer_n_5, "mu_nhit_layer_n_5");
    f.WriteTObject(&mu_nhit_layer_n_6, "mu_nhit_layer_n_6");
    f.WriteTObject(&mu_nhit_layer_n_7, "mu_nhit_layer_n_7");
    f.WriteTObject(&mu_nhit_layer_n_8, "mu_nhit_layer_n_8");
    f.WriteTObject(&mu_nhit_layer_n_9, "mu_nhit_layer_n_9");
    f.WriteTObject(&mu_nhit_layer_n_10, "mu_nhit_layer_n_10");
    f.WriteTObject(&mu_nhit_layer_n_11, "mu_nhit_layer_n_11");
    f.WriteTObject(&mu_nhit_layer_n_12, "mu_nhit_layer_n_12");
    f.WriteTObject(&mu_nhit_layer_n_13, "mu_nhit_layer_n_13");
    f.WriteTObject(&mu_nhit_layer_n_14, "mu_nhit_layer_n_14");
    f.WriteTObject(&sig_nhit_layer_n_0, "sig_nhit_layer_n_0");
    f.WriteTObject(&sig_nhit_layer_n_1, "sig_nhit_layer_n_1");
    f.WriteTObject(&sig_nhit_layer_n_2, "sig_nhit_layer_n_2");
    f.WriteTObject(&sig_nhit_layer_n_3, "sig_nhit_layer_n_3");
    f.WriteTObject(&sig_nhit_layer_n_4, "sig_nhit_layer_n_4");
    f.WriteTObject(&sig_nhit_layer_n_5, "sig_nhit_layer_n_5");
    f.WriteTObject(&sig_nhit_layer_n_6, "sig_nhit_layer_n_6");
    f.WriteTObject(&sig_nhit_layer_n_7, "sig_nhit_layer_n_7");
    f.WriteTObject(&sig_nhit_layer_n_8, "sig_nhit_layer_n_8");
    f.WriteTObject(&sig_nhit_layer_n_9, "sig_nhit_layer_n_9");
    f.WriteTObject(&sig_nhit_layer_n_10, "sig_nhit_layer_n_10");
    f.WriteTObject(&sig_nhit_layer_n_11, "sig_nhit_layer_n_11");
    f.WriteTObject(&sig_nhit_layer_n_12, "sig_nhit_layer_n_12");
    f.WriteTObject(&sig_nhit_layer_n_13, "sig_nhit_layer_n_13");
    f.WriteTObject(&sig_nhit_layer_n_14, "sig_nhit_layer_n_14");
    
    f.WriteTObject(&shower_nhit_max_layer, "shower_nhit_max_layer");
    f.WriteTObject(&shower_nhit_start_layer, "shower_nhit_start_layer");
    f.WriteTObject(&shower_nhit_end_layer, "shower_nhit_end_layer");
    f.WriteTObject(&shower_nhit_start_10_layer, "shower_nhit_start_10_layer");
    f.WriteTObject(&shower_nhit_end_10_layer, "shower_nhit_end_10_layer");
    f.WriteTObject(&shower_nhit_average, "shower_nhit_average");
    f.WriteTObject(&shower_nhit_max, "shower_nhit_max");
    
    f.WriteTObject(&mu_weight_layer_0, "mu_weight_layer_0");
    f.WriteTObject(&mu_weight_layer_1, "mu_weight_layer_1");
    f.WriteTObject(&mu_weight_layer_2, "mu_weight_layer_2");
    f.WriteTObject(&mu_weight_layer_3, "mu_weight_layer_3");
    f.WriteTObject(&mu_weight_layer_4, "mu_weight_layer_4");
    f.WriteTObject(&mu_weight_layer_5, "mu_weight_layer_5");
    f.WriteTObject(&mu_weight_layer_6, "mu_weight_layer_6");
    f.WriteTObject(&mu_weight_layer_7, "mu_weight_layer_7");
    f.WriteTObject(&mu_weight_layer_8, "mu_weight_layer_8");
    f.WriteTObject(&mu_weight_layer_9, "mu_weight_layer_9");
    f.WriteTObject(&mu_weight_layer_10, "mu_weight_layer_10");
    f.WriteTObject(&mu_weight_layer_11, "mu_weight_layer_11");
    f.WriteTObject(&mu_weight_layer_12, "mu_weight_layer_12");
    f.WriteTObject(&mu_weight_layer_13, "mu_weight_layer_13");
    f.WriteTObject(&mu_weight_layer_14, "mu_weight_layer_14");
    f.WriteTObject(&sig_weight_layer_0, "sig_weight_layer_0");
    f.WriteTObject(&sig_weight_layer_1, "sig_weight_layer_1");
    f.WriteTObject(&sig_weight_layer_2, "sig_weight_layer_2");
    f.WriteTObject(&sig_weight_layer_3, "sig_weight_layer_3");
    f.WriteTObject(&sig_weight_layer_4, "sig_weight_layer_4");
    f.WriteTObject(&sig_weight_layer_5, "sig_weight_layer_5");
    f.WriteTObject(&sig_weight_layer_6, "sig_weight_layer_6");
    f.WriteTObject(&sig_weight_layer_7, "sig_weight_layer_7");
    f.WriteTObject(&sig_weight_layer_8, "sig_weight_layer_8");
    f.WriteTObject(&sig_weight_layer_9, "sig_weight_layer_9");
    f.WriteTObject(&sig_weight_layer_10, "sig_weight_layer_10");
    f.WriteTObject(&sig_weight_layer_11, "sig_weight_layer_11");
    f.WriteTObject(&sig_weight_layer_12, "sig_weight_layer_12");
    f.WriteTObject(&sig_weight_layer_13, "sig_weight_layer_13");
    f.WriteTObject(&sig_weight_layer_14, "sig_weight_layer_14");

    f.WriteTObject(&mu_weight_layer_n_0, "mu_weight_layer_n_0");
    f.WriteTObject(&mu_weight_layer_n_1, "mu_weight_layer_n_1");
    f.WriteTObject(&mu_weight_layer_n_2, "mu_weight_layer_n_2");
    f.WriteTObject(&mu_weight_layer_n_3, "mu_weight_layer_n_3");
    f.WriteTObject(&mu_weight_layer_n_4, "mu_weight_layer_n_4");
    f.WriteTObject(&mu_weight_layer_n_5, "mu_weight_layer_n_5");
    f.WriteTObject(&mu_weight_layer_n_6, "mu_weight_layer_n_6");
    f.WriteTObject(&mu_weight_layer_n_7, "mu_weight_layer_n_7");
    f.WriteTObject(&mu_weight_layer_n_8, "mu_weight_layer_n_8");
    f.WriteTObject(&mu_weight_layer_n_9, "mu_weight_layer_n_9");
    f.WriteTObject(&mu_weight_layer_n_10, "mu_weight_layer_n_10");
    f.WriteTObject(&mu_weight_layer_n_11, "mu_weight_layer_n_11");
    f.WriteTObject(&mu_weight_layer_n_12, "mu_weight_layer_n_12");
    f.WriteTObject(&mu_weight_layer_n_13, "mu_weight_layer_n_13");
    f.WriteTObject(&mu_weight_layer_n_14, "mu_weight_layer_n_14");
    f.WriteTObject(&sig_weight_layer_n_0, "sig_weight_layer_n_0");
    f.WriteTObject(&sig_weight_layer_n_1, "sig_weight_layer_n_1");
    f.WriteTObject(&sig_weight_layer_n_2, "sig_weight_layer_n_2");
    f.WriteTObject(&sig_weight_layer_n_3, "sig_weight_layer_n_3");
    f.WriteTObject(&sig_weight_layer_n_4, "sig_weight_layer_n_4");
    f.WriteTObject(&sig_weight_layer_n_5, "sig_weight_layer_n_5");
    f.WriteTObject(&sig_weight_layer_n_6, "sig_weight_layer_n_6");
    f.WriteTObject(&sig_weight_layer_n_7, "sig_weight_layer_n_7");
    f.WriteTObject(&sig_weight_layer_n_8, "sig_weight_layer_n_8");
    f.WriteTObject(&sig_weight_layer_n_9, "sig_weight_layer_n_9");
    f.WriteTObject(&sig_weight_layer_n_10, "sig_weight_layer_n_10");
    f.WriteTObject(&sig_weight_layer_n_11, "sig_weight_layer_n_11");
    f.WriteTObject(&sig_weight_layer_n_12, "sig_weight_layer_n_12");
    f.WriteTObject(&sig_weight_layer_n_13, "sig_weight_layer_n_13");
    f.WriteTObject(&sig_weight_layer_n_14, "sig_weight_layer_n_14");
    
    f.WriteTObject(&shower_weight_max_layer, "shower_weight_max_layer");
    f.WriteTObject(&shower_weight_start_layer, "shower_weight_start_layer");
    f.WriteTObject(&shower_weight_end_layer, "shower_weight_end_layer");
    f.WriteTObject(&shower_weight_start_10_layer, "shower_weight_start_10_layer");
    f.WriteTObject(&shower_weight_end_10_layer, "shower_weight_end_10_layer");
    f.WriteTObject(&shower_weight_average, "shower_weight_average");
    f.WriteTObject(&shower_weight_max, "shower_weight_max");

    f.WriteTObject(&mu_barycenter_x_layer_0, "mu_barycenter_x_layer_0");
    f.WriteTObject(&mu_barycenter_x_layer_1, "mu_barycenter_x_layer_1");
    f.WriteTObject(&mu_barycenter_x_layer_2, "mu_barycenter_x_layer_2");
    f.WriteTObject(&mu_barycenter_x_layer_3, "mu_barycenter_x_layer_3");
    f.WriteTObject(&mu_barycenter_x_layer_4, "mu_barycenter_x_layer_4");
    f.WriteTObject(&mu_barycenter_x_layer_5, "mu_barycenter_x_layer_5");
    f.WriteTObject(&mu_barycenter_x_layer_6, "mu_barycenter_x_layer_6");
    f.WriteTObject(&mu_barycenter_x_layer_7, "mu_barycenter_x_layer_7");
    f.WriteTObject(&mu_barycenter_x_layer_8, "mu_barycenter_x_layer_8");
    f.WriteTObject(&mu_barycenter_x_layer_9, "mu_barycenter_x_layer_9");
    f.WriteTObject(&mu_barycenter_x_layer_10, "mu_barycenter_x_layer_10");
    f.WriteTObject(&mu_barycenter_x_layer_11, "mu_barycenter_x_layer_11");
    f.WriteTObject(&mu_barycenter_x_layer_12, "mu_barycenter_x_layer_12");
    f.WriteTObject(&mu_barycenter_x_layer_13, "mu_barycenter_x_layer_13");
    f.WriteTObject(&mu_barycenter_x_layer_14, "mu_barycenter_x_layer_14");
    f.WriteTObject(&sig_barycenter_x_layer_0, "sig_barycenter_x_layer_0");
    f.WriteTObject(&sig_barycenter_x_layer_1, "sig_barycenter_x_layer_1");
    f.WriteTObject(&sig_barycenter_x_layer_2, "sig_barycenter_x_layer_2");
    f.WriteTObject(&sig_barycenter_x_layer_3, "sig_barycenter_x_layer_3");
    f.WriteTObject(&sig_barycenter_x_layer_4, "sig_barycenter_x_layer_4");
    f.WriteTObject(&sig_barycenter_x_layer_5, "sig_barycenter_x_layer_5");
    f.WriteTObject(&sig_barycenter_x_layer_6, "sig_barycenter_x_layer_6");
    f.WriteTObject(&sig_barycenter_x_layer_7, "sig_barycenter_x_layer_7");
    f.WriteTObject(&sig_barycenter_x_layer_8, "sig_barycenter_x_layer_8");
    f.WriteTObject(&sig_barycenter_x_layer_9, "sig_barycenter_x_layer_9");
    f.WriteTObject(&sig_barycenter_x_layer_10, "sig_barycenter_x_layer_10");
    f.WriteTObject(&sig_barycenter_x_layer_11, "sig_barycenter_x_layer_11");
    f.WriteTObject(&sig_barycenter_x_layer_12, "sig_barycenter_x_layer_12");
    f.WriteTObject(&sig_barycenter_x_layer_13, "sig_barycenter_x_layer_13");
    f.WriteTObject(&sig_barycenter_x_layer_14, "sig_barycenter_x_layer_14");

    f.WriteTObject(&mu_barycenter_y_layer_0, "mu_barycenter_y_layer_0");
    f.WriteTObject(&mu_barycenter_y_layer_1, "mu_barycenter_y_layer_1");
    f.WriteTObject(&mu_barycenter_y_layer_2, "mu_barycenter_y_layer_2");
    f.WriteTObject(&mu_barycenter_y_layer_3, "mu_barycenter_y_layer_3");
    f.WriteTObject(&mu_barycenter_y_layer_4, "mu_barycenter_y_layer_4");
    f.WriteTObject(&mu_barycenter_y_layer_5, "mu_barycenter_y_layer_5");
    f.WriteTObject(&mu_barycenter_y_layer_6, "mu_barycenter_y_layer_6");
    f.WriteTObject(&mu_barycenter_y_layer_7, "mu_barycenter_y_layer_7");
    f.WriteTObject(&mu_barycenter_y_layer_8, "mu_barycenter_y_layer_8");
    f.WriteTObject(&mu_barycenter_y_layer_9, "mu_barycenter_y_layer_9");
    f.WriteTObject(&mu_barycenter_y_layer_10, "mu_barycenter_y_layer_10");
    f.WriteTObject(&mu_barycenter_y_layer_11, "mu_barycenter_y_layer_11");
    f.WriteTObject(&mu_barycenter_y_layer_12, "mu_barycenter_y_layer_12");
    f.WriteTObject(&mu_barycenter_y_layer_13, "mu_barycenter_y_layer_13");
    f.WriteTObject(&mu_barycenter_y_layer_14, "mu_barycenter_y_layer_14");
    f.WriteTObject(&sig_barycenter_y_layer_0, "sig_barycenter_y_layer_0");
    f.WriteTObject(&sig_barycenter_y_layer_1, "sig_barycenter_y_layer_1");
    f.WriteTObject(&sig_barycenter_y_layer_2, "sig_barycenter_y_layer_2");
    f.WriteTObject(&sig_barycenter_y_layer_3, "sig_barycenter_y_layer_3");
    f.WriteTObject(&sig_barycenter_y_layer_4, "sig_barycenter_y_layer_4");
    f.WriteTObject(&sig_barycenter_y_layer_5, "sig_barycenter_y_layer_5");
    f.WriteTObject(&sig_barycenter_y_layer_6, "sig_barycenter_y_layer_6");
    f.WriteTObject(&sig_barycenter_y_layer_7, "sig_barycenter_y_layer_7");
    f.WriteTObject(&sig_barycenter_y_layer_8, "sig_barycenter_y_layer_8");
    f.WriteTObject(&sig_barycenter_y_layer_9, "sig_barycenter_y_layer_9");
    f.WriteTObject(&sig_barycenter_y_layer_10, "sig_barycenter_y_layer_10");
    f.WriteTObject(&sig_barycenter_y_layer_11, "sig_barycenter_y_layer_11");
    f.WriteTObject(&sig_barycenter_y_layer_12, "sig_barycenter_y_layer_12");
    f.WriteTObject(&sig_barycenter_y_layer_13, "sig_barycenter_y_layer_13");
    f.WriteTObject(&sig_barycenter_y_layer_14, "sig_barycenter_y_layer_14");
    f.WriteTObject(&mu_barycenter_r_layer_0, "mu_barycenter_r_layer_0");
    f.WriteTObject(&mu_barycenter_r_layer_1, "mu_barycenter_r_layer_1");
    f.WriteTObject(&mu_barycenter_r_layer_2, "mu_barycenter_r_layer_2");
    f.WriteTObject(&mu_barycenter_r_layer_3, "mu_barycenter_r_layer_3");
    f.WriteTObject(&mu_barycenter_r_layer_4, "mu_barycenter_r_layer_4");
    f.WriteTObject(&mu_barycenter_r_layer_5, "mu_barycenter_r_layer_5");
    f.WriteTObject(&mu_barycenter_r_layer_6, "mu_barycenter_r_layer_6");
    f.WriteTObject(&mu_barycenter_r_layer_7, "mu_barycenter_r_layer_7");
    f.WriteTObject(&mu_barycenter_r_layer_8, "mu_barycenter_r_layer_8");
    f.WriteTObject(&mu_barycenter_r_layer_9, "mu_barycenter_r_layer_9");
    f.WriteTObject(&mu_barycenter_r_layer_10, "mu_barycenter_r_layer_10");
    f.WriteTObject(&mu_barycenter_r_layer_11, "mu_barycenter_r_layer_11");
    f.WriteTObject(&mu_barycenter_r_layer_12, "mu_barycenter_r_layer_12");
    f.WriteTObject(&mu_barycenter_r_layer_13, "mu_barycenter_r_layer_13");
    f.WriteTObject(&mu_barycenter_r_layer_14, "mu_barycenter_r_layer_14");
    f.WriteTObject(&sig_barycenter_r_layer_0, "sig_barycenter_r_layer_0");
    f.WriteTObject(&sig_barycenter_r_layer_1, "sig_barycenter_r_layer_1");
    f.WriteTObject(&sig_barycenter_r_layer_2, "sig_barycenter_r_layer_2");
    f.WriteTObject(&sig_barycenter_r_layer_3, "sig_barycenter_r_layer_3");
    f.WriteTObject(&sig_barycenter_r_layer_4, "sig_barycenter_r_layer_4");
    f.WriteTObject(&sig_barycenter_r_layer_5, "sig_barycenter_r_layer_5");
    f.WriteTObject(&sig_barycenter_r_layer_6, "sig_barycenter_r_layer_6");
    f.WriteTObject(&sig_barycenter_r_layer_7, "sig_barycenter_r_layer_7");
    f.WriteTObject(&sig_barycenter_r_layer_8, "sig_barycenter_r_layer_8");
    f.WriteTObject(&sig_barycenter_r_layer_9, "sig_barycenter_r_layer_9");
    f.WriteTObject(&sig_barycenter_r_layer_10, "sig_barycenter_r_layer_10");
    f.WriteTObject(&sig_barycenter_r_layer_11, "sig_barycenter_r_layer_11");
    f.WriteTObject(&sig_barycenter_r_layer_12, "sig_barycenter_r_layer_12");
    f.WriteTObject(&sig_barycenter_r_layer_13, "sig_barycenter_r_layer_13");
    f.WriteTObject(&sig_barycenter_r_layer_14, "sig_barycenter_r_layer_14");

    f.WriteTObject(&mu_barycenter_x, "mu_barycenter_x");
    f.WriteTObject(&mu_barycenter_y, "mu_barycenter_y");
    f.WriteTObject(&mu_barycenter_z, "mu_barycenter_z");
    f.WriteTObject(&mu_barycenter_r, "mu_barycenter_r");
    f.WriteTObject(&sig_barycenter_x, "sig_barycenter_x");
    f.WriteTObject(&sig_barycenter_y, "sig_barycenter_y");
    f.WriteTObject(&sig_barycenter_z, "sig_barycenter_z");
    f.WriteTObject(&sig_barycenter_r, "sig_barycenter_r");

    f.WriteObject(&energies, "energies");
    f.WriteObject(&energies_tr, "energies_tr");
    f.WriteObject(&W_thicknesses, "W_thicknesses");
    f.Close();

}
