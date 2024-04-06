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

#define N_ENERGIES 3 //18
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

void graph_setup_add(TGraph *g, string title, Color_t color){
    g->SetTitle(title.c_str());
    g->SetLineColor(color);
    g->SetLineWidth(3);
    return;
}


void analysis (string particle, bool masked=false) {
    
    // double energies[N_ENERGIES] = {1, 2, 5, 8, 10, 20, 40, 60, 80, 100, 150};
    double test_e[N_ENERGIES]={6., 60., 150};
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
    
    TString result_name = "resolution_"+particle+"_result.root" ;
    
    TFile f(result_name, "recreate");
    TTree *outtree = new TTree("ntp","NTuples");

    // branches definitions
    // We will save both the branches and the histos
    // Values
    Float_t b_nhit, b_sume, b_weighte, b_bar_x, b_bar_y, b_bar_z, b_bar_r;
    Float_t b_mol;
    Float_t b_bar_x_layer_0, b_bar_x_layer_1, b_bar_x_layer_2, b_bar_x_layer_3, b_bar_x_layer_4, b_bar_x_layer_5, b_bar_x_layer_6, b_bar_x_layer_7, b_bar_x_layer_8, b_bar_x_layer_9, b_bar_x_layer_10, b_bar_x_layer_11, b_bar_x_layer_12, b_bar_x_layer_13, b_bar_x_layer_14;
    Float_t b_bar_y_layer_0, b_bar_y_layer_1, b_bar_y_layer_2, b_bar_y_layer_3, b_bar_y_layer_4, b_bar_y_layer_5, b_bar_y_layer_6, b_bar_y_layer_7, b_bar_y_layer_8, b_bar_y_layer_9, b_bar_y_layer_10, b_bar_y_layer_11, b_bar_y_layer_12, b_bar_y_layer_13, b_bar_y_layer_14;
    Float_t b_bar_r_layer_0, b_bar_r_layer_1, b_bar_r_layer_2, b_bar_r_layer_3, b_bar_r_layer_4, b_bar_r_layer_5, b_bar_r_layer_6, b_bar_r_layer_7, b_bar_r_layer_8, b_bar_r_layer_9, b_bar_r_layer_10, b_bar_r_layer_11, b_bar_r_layer_12, b_bar_r_layer_13, b_bar_r_layer_14;
    Float_t b_shower_nhit_max_layer, b_shower_nhit_start_layer, b_shower_nhit_end_layer, b_shower_nhit_start_10_layer, b_shower_nhit_end_10_layer, b_shower_nhit_average, b_shower_nhit_max;
    Float_t b_shower_weighte_max_layer, b_shower_weighte_start_layer, b_shower_weighte_end_layer, b_shower_weighte_start_10_layer, b_shower_weighte_end_10_layer, b_shower_weighte_average, b_shower_weighte_max;
    Float_t b_nhit_layer_0, b_nhit_layer_1, b_nhit_layer_2, b_nhit_layer_3, b_nhit_layer_4, b_nhit_layer_5, b_nhit_layer_6, b_nhit_layer_7, b_nhit_layer_8, b_nhit_layer_9, b_nhit_layer_10, b_nhit_layer_11, b_nhit_layer_12, b_nhit_layer_13, b_nhit_layer_14;
    Float_t b_nhit_layer_n_0, b_nhit_layer_n_1, b_nhit_layer_n_2, b_nhit_layer_n_3, b_nhit_layer_n_4, b_nhit_layer_n_5, b_nhit_layer_n_6, b_nhit_layer_n_7, b_nhit_layer_n_8, b_nhit_layer_n_9, b_nhit_layer_n_10, b_nhit_layer_n_11, b_nhit_layer_n_12, b_nhit_layer_n_13, b_nhit_layer_n_14;
    Float_t b_weighte_layer_0, b_weighte_layer_1, b_weighte_layer_2, b_weighte_layer_3, b_weighte_layer_4, b_weighte_layer_5, b_weighte_layer_6, b_weighte_layer_7, b_weighte_layer_8, b_weighte_layer_9, b_weighte_layer_10, b_weighte_layer_11, b_weighte_layer_12, b_weighte_layer_13, b_weighte_layer_14;
    Float_t b_weighte_layer_n_0, b_weighte_layer_n_1, b_weighte_layer_n_2, b_weighte_layer_n_3, b_weighte_layer_n_4, b_weighte_layer_n_5, b_weighte_layer_n_6, b_weighte_layer_n_7, b_weighte_layer_n_8, b_weighte_layer_n_9, b_weighte_layer_n_10, b_weighte_layer_n_11, b_weighte_layer_n_12, b_weighte_layer_n_13, b_weighte_layer_n_14;

    // Adding the branches
    outtree->Branch("nhit",&b_nhit,"b_nhit/F");
    outtree->Branch("sume",&b_sume,"b_sume/F");
    outtree->Branch("weighte",&b_weighte,"b_weighte/F");
    outtree->Branch("mol",&b_mol,"mol/F");
    outtree->Branch("bar_x",&b_bar_x,"b_bar_x/F");
    outtree->Branch("bar_y",&b_bar_y,"b_bar_y/F");
    outtree->Branch("bar_z",&b_bar_z,"b_bar_z/F");
    outtree->Branch("bar_r",&b_bar_r,"b_bar_r/F");
    outtree->Branch("bar_x_layer_0",&b_bar_x_layer_0,"b_bar_x_layer_0/F");
    outtree->Branch("bar_x_layer_1",&b_bar_x_layer_1,"b_bar_x_layer_1/F");
    outtree->Branch("bar_x_layer_2",&b_bar_x_layer_2,"b_bar_x_layer_2/F");
    outtree->Branch("bar_x_layer_3",&b_bar_x_layer_3,"b_bar_x_layer_3/F");
    outtree->Branch("bar_x_layer_4",&b_bar_x_layer_4,"b_bar_x_layer_4/F");
    outtree->Branch("bar_x_layer_5",&b_bar_x_layer_5,"b_bar_x_layer_5/F");
    outtree->Branch("bar_x_layer_6",&b_bar_x_layer_6,"b_bar_x_layer_6/F");
    outtree->Branch("bar_x_layer_7",&b_bar_x_layer_7,"b_bar_x_layer_7/F");
    outtree->Branch("bar_x_layer_8",&b_bar_x_layer_8,"b_bar_x_layer_8/F");
    outtree->Branch("bar_x_layer_9",&b_bar_x_layer_9,"b_bar_x_layer_9/F");
    outtree->Branch("bar_x_layer_10",&b_bar_x_layer_10,"b_bar_x_layer_10/F");
    outtree->Branch("bar_x_layer_11",&b_bar_x_layer_11,"b_bar_x_layer_11/F");
    outtree->Branch("bar_x_layer_12",&b_bar_x_layer_12,"b_bar_x_layer_12/F");
    outtree->Branch("bar_x_layer_13",&b_bar_x_layer_13,"b_bar_x_layer_13/F");
    outtree->Branch("bar_x_layer_14",&b_bar_x_layer_14,"b_bar_x_layer_14/F");
    outtree->Branch("bar_y_layer_0",&b_bar_y_layer_0,"b_bar_y_layer_0/F");
    outtree->Branch("bar_y_layer_1",&b_bar_y_layer_1,"b_bar_y_layer_1/F");
    outtree->Branch("bar_y_layer_2",&b_bar_y_layer_2,"b_bar_y_layer_2/F");
    outtree->Branch("bar_y_layer_3",&b_bar_y_layer_3,"b_bar_y_layer_3/F");
    outtree->Branch("bar_y_layer_4",&b_bar_y_layer_4,"b_bar_y_layer_4/F");
    outtree->Branch("bar_y_layer_5",&b_bar_y_layer_5,"b_bar_y_layer_5/F");
    outtree->Branch("bar_y_layer_6",&b_bar_y_layer_6,"b_bar_y_layer_6/F");
    outtree->Branch("bar_y_layer_7",&b_bar_y_layer_7,"b_bar_y_layer_7/F");
    outtree->Branch("bar_y_layer_8",&b_bar_y_layer_8,"b_bar_y_layer_8/F");
    outtree->Branch("bar_y_layer_9",&b_bar_y_layer_9,"b_bar_y_layer_9/F");
    outtree->Branch("bar_y_layer_10",&b_bar_y_layer_10,"b_bar_y_layer_10/F");
    outtree->Branch("bar_y_layer_11",&b_bar_y_layer_11,"b_bar_y_layer_11/F");
    outtree->Branch("bar_y_layer_12",&b_bar_y_layer_12,"b_bar_y_layer_12/F");
    outtree->Branch("bar_y_layer_13",&b_bar_y_layer_13,"b_bar_y_layer_13/F");
    outtree->Branch("bar_y_layer_14",&b_bar_y_layer_14,"b_bar_y_layer_14/F");
    outtree->Branch("bar_r_layer_0",&b_bar_r_layer_0,"b_bar_r_layer_0/F");
    outtree->Branch("bar_r_layer_1",&b_bar_r_layer_1,"b_bar_r_layer_1/F");
    outtree->Branch("bar_r_layer_2",&b_bar_r_layer_2,"b_bar_r_layer_2/F");
    outtree->Branch("bar_r_layer_3",&b_bar_r_layer_3,"b_bar_r_layer_3/F");
    outtree->Branch("bar_r_layer_4",&b_bar_r_layer_4,"b_bar_r_layer_4/F");
    outtree->Branch("bar_r_layer_5",&b_bar_r_layer_5,"b_bar_r_layer_5/F");
    outtree->Branch("bar_r_layer_6",&b_bar_r_layer_6,"b_bar_r_layer_6/F");
    outtree->Branch("bar_r_layer_7",&b_bar_r_layer_7,"b_bar_r_layer_7/F");
    outtree->Branch("bar_r_layer_8",&b_bar_r_layer_8,"b_bar_r_layer_8/F");
    outtree->Branch("bar_r_layer_9",&b_bar_r_layer_9,"b_bar_r_layer_9/F");
    outtree->Branch("bar_r_layer_10",&b_bar_r_layer_10,"b_bar_r_layer_10/F");
    outtree->Branch("bar_r_layer_11",&b_bar_r_layer_11,"b_bar_r_layer_11/F");
    outtree->Branch("bar_r_layer_12",&b_bar_r_layer_12,"b_bar_r_layer_12/F");
    outtree->Branch("bar_r_layer_13",&b_bar_r_layer_13,"b_bar_r_layer_13/F");
    outtree->Branch("bar_r_layer_14",&b_bar_r_layer_14,"b_bar_r_layer_14/F");
    outtree->Branch("shower_nhit_max_layer",&b_shower_nhit_max_layer,"b_shower_nhit_max_layer/F");
    outtree->Branch("shower_nhit_start_layer",&b_shower_nhit_start_layer,"b_shower_nhit_start_layer/F");
    outtree->Branch("shower_nhit_end_layer",&b_shower_nhit_end_layer,"b_shower_nhit_end_layer/F");
    outtree->Branch("shower_nhit_start_10_layer",&b_shower_nhit_start_10_layer,"b_shower_nhit_start_10_layer/F");
    outtree->Branch("shower_nhit_end_10_layer",&b_shower_nhit_end_10_layer,"b_shower_nhit_end_10_layer/F");
    outtree->Branch("shower_nhit_average",&b_shower_nhit_average,"b_shower_nhit_average/F");
    outtree->Branch("shower_nhit_max",&b_shower_nhit_max,"b_shower_nhit_max/F");
    outtree->Branch("shower_weighte_max_layer",&b_shower_weighte_max_layer,"b_shower_weighte_max_layer/F");
    outtree->Branch("shower_weighte_start_layer",&b_shower_weighte_start_layer,"b_shower_weighte_start_layer/F");
    outtree->Branch("shower_weighte_end_layer",&b_shower_weighte_end_layer,"b_shower_weighte_end_layer/F");
    outtree->Branch("shower_weighte_start_10_layer",&b_shower_weighte_start_10_layer,"b_shower_weighte_start_10_layer/F");
    outtree->Branch("shower_weighte_end_10_layer",&b_shower_weighte_end_10_layer,"b_shower_weighte_end_10_layer/F");
    outtree->Branch("shower_weighte_average",&b_shower_weighte_average,"b_shower_weighte_average/F");
    outtree->Branch("shower_weighte_max",&b_shower_weighte_max,"b_shower_weighte_max/F");
    outtree->Branch("nhit_layer_0",&b_nhit_layer_0,"b_nhit_layer_0/F");
    outtree->Branch("nhit_layer_1",&b_nhit_layer_1,"b_nhit_layer_1/F");
    outtree->Branch("nhit_layer_2",&b_nhit_layer_2,"b_nhit_layer_2/F");
    outtree->Branch("nhit_layer_3",&b_nhit_layer_3,"b_nhit_layer_3/F");
    outtree->Branch("nhit_layer_4",&b_nhit_layer_4,"b_nhit_layer_4/F");
    outtree->Branch("nhit_layer_5",&b_nhit_layer_5,"b_nhit_layer_5/F");
    outtree->Branch("nhit_layer_6",&b_nhit_layer_6,"b_nhit_layer_6/F");
    outtree->Branch("nhit_layer_7",&b_nhit_layer_7,"b_nhit_layer_7/F");
    outtree->Branch("nhit_layer_8",&b_nhit_layer_8,"b_nhit_layer_8/F");
    outtree->Branch("nhit_layer_9",&b_nhit_layer_9,"b_nhit_layer_9/F");
    outtree->Branch("nhit_layer_10",&b_nhit_layer_10,"b_nhit_layer_10/F");
    outtree->Branch("nhit_layer_11",&b_nhit_layer_11,"b_nhit_layer_11/F");
    outtree->Branch("nhit_layer_12",&b_nhit_layer_12,"b_nhit_layer_12/F");
    outtree->Branch("nhit_layer_13",&b_nhit_layer_13,"b_nhit_layer_13/F");
    outtree->Branch("nhit_layer_14",&b_nhit_layer_14,"b_nhit_layer_14/F");
    outtree->Branch("nhit_layer_n_0",&b_nhit_layer_n_0,"b_nhit_layer_n_0/F");
    outtree->Branch("nhit_layer_n_1",&b_nhit_layer_n_1,"b_nhit_layer_n_1/F");
    outtree->Branch("nhit_layer_n_2",&b_nhit_layer_n_2,"b_nhit_layer_n_2/F");
    outtree->Branch("nhit_layer_n_3",&b_nhit_layer_n_3,"b_nhit_layer_n_3/F");
    outtree->Branch("nhit_layer_n_4",&b_nhit_layer_n_4,"b_nhit_layer_n_4/F");
    outtree->Branch("nhit_layer_n_5",&b_nhit_layer_n_5,"b_nhit_layer_n_5/F");
    outtree->Branch("nhit_layer_n_6",&b_nhit_layer_n_6,"b_nhit_layer_n_6/F");
    outtree->Branch("nhit_layer_n_7",&b_nhit_layer_n_7,"b_nhit_layer_n_7/F");
    outtree->Branch("nhit_layer_n_8",&b_nhit_layer_n_8,"b_nhit_layer_n_8/F");
    outtree->Branch("nhit_layer_n_9",&b_nhit_layer_n_9,"b_nhit_layer_n_9/F");
    outtree->Branch("nhit_layer_n_10",&b_nhit_layer_n_10,"b_nhit_layer_n_10/F");
    outtree->Branch("nhit_layer_n_11",&b_nhit_layer_n_11,"b_nhit_layer_n_11/F");
    outtree->Branch("nhit_layer_n_12",&b_nhit_layer_n_12,"b_nhit_layer_n_12/F");
    outtree->Branch("nhit_layer_n_13",&b_nhit_layer_n_13,"b_nhit_layer_n_13/F");
    outtree->Branch("nhit_layer_n_14",&b_nhit_layer_n_14,"b_nhit_layer_n_14/F");
    outtree->Branch("weighte_layer_0",&b_weighte_layer_0,"b_weighte_layer_0/F");
    outtree->Branch("weighte_layer_1",&b_weighte_layer_1,"b_weighte_layer_1/F");
    outtree->Branch("weighte_layer_2",&b_weighte_layer_2,"b_weighte_layer_2/F");
    outtree->Branch("weighte_layer_3",&b_weighte_layer_3,"b_weighte_layer_3/F");
    outtree->Branch("weighte_layer_4",&b_weighte_layer_4,"b_weighte_layer_4/F");
    outtree->Branch("weighte_layer_5",&b_weighte_layer_5,"b_weighte_layer_5/F");
    outtree->Branch("weighte_layer_6",&b_weighte_layer_6,"b_weighte_layer_6/F");
    outtree->Branch("weighte_layer_7",&b_weighte_layer_7,"b_weighte_layer_7/F");
    outtree->Branch("weighte_layer_8",&b_weighte_layer_8,"b_weighte_layer_8/F");
    outtree->Branch("weighte_layer_9",&b_weighte_layer_9,"b_weighte_layer_9/F");
    outtree->Branch("weighte_layer_10",&b_weighte_layer_10,"b_weighte_layer_10/F");
    outtree->Branch("weighte_layer_11",&b_weighte_layer_11,"b_weighte_layer_11/F");
    outtree->Branch("weighte_layer_12",&b_weighte_layer_12,"b_weighte_layer_12/F");
    outtree->Branch("weighte_layer_13",&b_weighte_layer_13,"b_weighte_layer_13/F");
    outtree->Branch("weighte_layer_14",&b_weighte_layer_14,"b_weighte_layer_14/F");
    outtree->Branch("weighte_layer_n_0",&b_weighte_layer_n_0,"b_weighte_layer_n_0/F");
    outtree->Branch("weighte_layer_n_1",&b_weighte_layer_n_1,"b_weighte_layer_n_1/F");
    outtree->Branch("weighte_layer_n_2",&b_weighte_layer_n_2,"b_weighte_layer_n_2/F");
    outtree->Branch("weighte_layer_n_3",&b_weighte_layer_n_3,"b_weighte_layer_n_3/F");
    outtree->Branch("weighte_layer_n_4",&b_weighte_layer_n_4,"b_weighte_layer_n_4/F");
    outtree->Branch("weighte_layer_n_5",&b_weighte_layer_n_5,"b_weighte_layer_n_5/F");
    outtree->Branch("weighte_layer_n_6",&b_weighte_layer_n_6,"b_weighte_layer_n_6/F");
    outtree->Branch("weighte_layer_n_7",&b_weighte_layer_n_7,"b_weighte_layer_n_7/F");
    outtree->Branch("weighte_layer_n_8",&b_weighte_layer_n_8,"b_weighte_layer_n_8/F");
    outtree->Branch("weighte_layer_n_9",&b_weighte_layer_n_9,"b_weighte_layer_n_9/F");
    outtree->Branch("weighte_layer_n_10",&b_weighte_layer_n_10,"b_weighte_layer_n_10/F");
    outtree->Branch("weighte_layer_n_11",&b_weighte_layer_n_11,"b_weighte_layer_n_11/F");
    outtree->Branch("weighte_layer_n_12",&b_weighte_layer_n_12,"b_weighte_layer_n_12/F");
    outtree->Branch("weighte_layer_n_13",&b_weighte_layer_n_13,"b_weighte_layer_n_13/F");
    outtree->Branch("weighte_layer_n_14",&b_weighte_layer_n_14,"b_weighte_layer_n_14/F");


    // Histo definitions, before the energy loop
    
    // Resolution histos
    string part_string = particle;

    TH1F *h_nhit = new TH1F(("NumHits_" + part_string ).c_str(), ("Number of hits " + part_string ).c_str(), 4001, -0.5, 4000.5);
    TH1F *h_sume = new TH1F(("SumEnergy_" + part_string ).c_str(), ("Sum Energy " + part_string ).c_str(), 5001, -0.5, 15000.5);
    TH1F *h_weighte = new TH1F(("WSumEnergy_" + part_string ).c_str(), ("W Sum Energy " + part_string ).c_str(), 5001, -0.5, 50000.5);

    // Barycenter histos
    TH1F *h_bar_x = new TH1F(("Barycenter_x_" + part_string ).c_str(), ("Barycenter (x-axis) " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y = new TH1F(("Barycenter_y_" + part_string ).c_str(), ("Barycenter (y-axis) " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_z = new TH1F(("Barycenter_z_" + part_string ).c_str(), ("Barycenter (z-axis) " + part_string ).c_str(), 1000, -500., 500.);
    TH1F *h_bar_r = new TH1F(("Barycenter_r_" + part_string ).c_str(), ("Barycenter (radial) " + part_string ).c_str(), 201, -100.5, 100.5);

    // Barycenter per layer
    TH1F *h_bar_r_layer_0 = new TH1F(("Barycenter_r_layer_0_" + part_string ).c_str(), ("Barycenter (radial) layer 0 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_1 = new TH1F(("Barycenter_r_layer_1_" + part_string ).c_str(), ("Barycenter (radial) layer 1 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_2 = new TH1F(("Barycenter_r_layer_2_" + part_string ).c_str(), ("Barycenter (radial) layer 2 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_3 = new TH1F(("Barycenter_r_layer_3_" + part_string ).c_str(), ("Barycenter (radial) layer 3 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_4 = new TH1F(("Barycenter_r_layer_4_" + part_string ).c_str(), ("Barycenter (radial) layer 4 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_5 = new TH1F(("Barycenter_r_layer_5_" + part_string ).c_str(), ("Barycenter (radial) layer 5 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_6 = new TH1F(("Barycenter_r_layer_6_" + part_string ).c_str(), ("Barycenter (radial) layer 6 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_7 = new TH1F(("Barycenter_r_layer_7_" + part_string ).c_str(), ("Barycenter (radial) layer 7 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_8 = new TH1F(("Barycenter_r_layer_8_" + part_string ).c_str(), ("Barycenter (radial) layer 8 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_9 = new TH1F(("Barycenter_r_layer_9_" + part_string ).c_str(), ("Barycenter (radial) layer 9 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_10 = new TH1F(("Barycenter_r_layer_10_" + part_string ).c_str(), ("Barycenter (radial) layer 10 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_11 = new TH1F(("Barycenter_r_layer_11_" + part_string ).c_str(), ("Barycenter (radial) layer 11 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_12 = new TH1F(("Barycenter_r_layer_12_" + part_string ).c_str(), ("Barycenter (radial) layer 12 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_13 = new TH1F(("Barycenter_r_layer_13_" + part_string ).c_str(), ("Barycenter (radial) layer 13 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_r_layer_14 = new TH1F(("Barycenter_r_layer_14_" + part_string ).c_str(), ("Barycenter (radial) layer 14 " + part_string ).c_str(), 201, -100.5, 100.5);

    TH1F *h_bar_x_layer_0 = new TH1F(("Barycenter_x_layer_0_" + part_string ).c_str(), ("Barycenter (x-axis) layer 0 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_1 = new TH1F(("Barycenter_x_layer_1_" + part_string ).c_str(), ("Barycenter (x-axis) layer 1 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_2 = new TH1F(("Barycenter_x_layer_2_" + part_string ).c_str(), ("Barycenter (x-axis) layer 2 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_3 = new TH1F(("Barycenter_x_layer_3_" + part_string ).c_str(), ("Barycenter (x-axis) layer 3 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_4 = new TH1F(("Barycenter_x_layer_4_" + part_string ).c_str(), ("Barycenter (x-axis) layer 4 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_5 = new TH1F(("Barycenter_x_layer_5_" + part_string ).c_str(), ("Barycenter (x-axis) layer 5 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_6 = new TH1F(("Barycenter_x_layer_6_" + part_string ).c_str(), ("Barycenter (x-axis) layer 6 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_7 = new TH1F(("Barycenter_x_layer_7_" + part_string ).c_str(), ("Barycenter (x-axis) layer 7 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_8 = new TH1F(("Barycenter_x_layer_8_" + part_string ).c_str(), ("Barycenter (x-axis) layer 8 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_9 = new TH1F(("Barycenter_x_layer_9_" + part_string ).c_str(), ("Barycenter (x-axis) layer 9 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_10 = new TH1F(("Barycenter_x_layer_10_" + part_string ).c_str(), ("Barycenter (x-axis) layer 10 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_11 = new TH1F(("Barycenter_x_layer_11_" + part_string ).c_str(), ("Barycenter (x-axis) layer 11 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_12 = new TH1F(("Barycenter_x_layer_12_" + part_string ).c_str(), ("Barycenter (x-axis) layer 12 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_13 = new TH1F(("Barycenter_x_layer_13_" + part_string ).c_str(), ("Barycenter (x-axis) layer 13 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_x_layer_14 = new TH1F(("Barycenter_x_layer_14_" + part_string ).c_str(), ("Barycenter (x-axis) layer 14 " + part_string ).c_str(), 201, -100.5, 100.5);

    TH1F *h_bar_y_layer_0 = new TH1F(("Barycenter_y_layer_0_" + part_string ).c_str(), ("Barycenter (y-axis) layer 0 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_1 = new TH1F(("Barycenter_y_layer_1_" + part_string ).c_str(), ("Barycenter (y-axis) layer 1 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_2 = new TH1F(("Barycenter_y_layer_2_" + part_string ).c_str(), ("Barycenter (x-axis) layer 2 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_3 = new TH1F(("Barycenter_y_layer_3_" + part_string ).c_str(), ("Barycenter (y-axis) layer 3 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_4 = new TH1F(("Barycenter_y_layer_4_" + part_string ).c_str(), ("Barycenter (y-axis) layer 4 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_5 = new TH1F(("Barycenter_y_layer_5_" + part_string ).c_str(), ("Barycenter (y-axis) layer 5 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_6 = new TH1F(("Barycenter_y_layer_6_" + part_string ).c_str(), ("Barycenter (y-axis) layer 6 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_7 = new TH1F(("Barycenter_y_layer_7_" + part_string ).c_str(), ("Barycenter (y-axis) layer 7 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_8 = new TH1F(("Barycenter_y_layer_8_" + part_string ).c_str(), ("Barycenter (y-axis) layer 8 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_9 = new TH1F(("Barycenter_y_layer_9_" + part_string ).c_str(), ("Barycenter (y-axis) layer 9 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_10 = new TH1F(("Barycenter_y_layer_10_" + part_string ).c_str(), ("Barycenter (y-axis) layer 10 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_11 = new TH1F(("Barycenter_y_layer_11_" + part_string ).c_str(), ("Barycenter (y-axis) layer 11 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_12 = new TH1F(("Barycenter_y_layer_12_" + part_string ).c_str(), ("Barycenter (y-axis) layer 12 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_13 = new TH1F(("Barycenter_y_layer_13_" + part_string ).c_str(), ("Barycenter (y-axis) layer 13 " + part_string ).c_str(), 201, -100.5, 100.5);
    TH1F *h_bar_y_layer_14 = new TH1F(("Barycenter_y_layer_14_" + part_string ).c_str(), ("Barycenter (y-axis) layer 14 " + part_string ).c_str(), 201, -100.5, 100.5);

    // Moliere histos
    TH1F *h_mol = new TH1F(("Radius90_" + part_string ).c_str(), ("Radius containing 90% of energy " + part_string ).c_str(), 101, -0.5, 100.5);

    // Shower profile characteristics
    TH1F *h_shower_nhit_max_layer = new TH1F(("ShowerNhitMaxLayer_" + part_string ).c_str(), ("Shower Nhit Max. (layer) " + part_string ).c_str(), 16, -1.5, 14.5);
    TH1F *h_shower_nhit_start_layer = new TH1F(("ShowerNhitStartLayer_" + part_string ).c_str(), ("Shower Nhit Start (layer) " + part_string ).c_str(), 16, -1.5, 14.5);
    TH1F *h_shower_nhit_end_layer = new TH1F(("ShowerNhitEndLayer_" + part_string ).c_str(), ("Shower Nhit End (layer) " + part_string ).c_str(), 16, -1.5, 14.5);
    TH1F *h_shower_nhit_start_10_layer = new TH1F(("ShowerNhitStart10Layer_" + part_string ).c_str(), ("Shower Nhit Start 0.1*Max (layer) " + part_string ).c_str(), 16, -1.5, 14.5);
    TH1F *h_shower_nhit_end_10_layer = new TH1F(("ShowerNhitEnd10Layer_" + part_string ).c_str(), ("Shower Nhit End 0.1*Max (layer) " + part_string ).c_str(), 16, -1.5, 14.5);
    TH1F *h_shower_nhit_average= new TH1F(("ShowerNhitAverage_" + part_string ).c_str(), ("Shower Nhit Average  " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_shower_nhit_max= new TH1F(("ShowerNhitMax_" + part_string ).c_str(), ("Shower Nhit Max  " + part_string ).c_str(), 201, -0.5, 200.5);

    TH1F *h_shower_weighte_max_layer = new TH1F(("ShowerWeightEMaxLayer_" + part_string ).c_str(), ("Shower WeightE Max. (layer) " + part_string ).c_str(), 16, -1.5, 14.5);
    TH1F *h_shower_weighte_start_layer = new TH1F(("ShowerWeightEStartLayer_" + part_string ).c_str(), ("Shower WeightE Start (layer) " + part_string ).c_str(), 16, -1.5, 14.5);
    TH1F *h_shower_weighte_end_layer = new TH1F(("ShowerWeightEEndLayer_" + part_string ).c_str(), ("Shower WeightE End (layer) " + part_string ).c_str(), 16, -1.5, 14.5);
    TH1F *h_shower_weighte_start_10_layer = new TH1F(("ShowerWeightEStart10Layer_" + part_string ).c_str(), ("Shower WeightE Start 0.1*Max (layer) " + part_string ).c_str(), 16, -1.5, 14.5);
    TH1F *h_shower_weighte_end_10_layer = new TH1F(("ShowerWeightEEnd10Layer_" + part_string ).c_str(), ("Shower WeightE End 0.1*Max (layer) " + part_string ).c_str(), 16, -1.5, 14.5);
    TH1F *h_shower_weighte_average= new TH1F(("ShowerWeightEAverage_" + part_string ).c_str(), ("Shower WeightE Average  " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_shower_weighte_max= new TH1F(("ShowerWeightEMax_" + part_string ).c_str(), ("Shower WeightE Max  " + part_string ).c_str(), 201, -0.5, 200.5);

    // Shower profile per layer
    TH1F *h_nhit_layer_0 = new TH1F(("NumHits_layer_0_" + part_string ).c_str(), ("Number of hits (layer 0) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_1 = new TH1F(("NumHits_layer_1_" + part_string ).c_str(), ("Number of hits (layer 1) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_2 = new TH1F(("NumHits_layer_2_" + part_string ).c_str(), ("Number of hits (layer 2) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_3 = new TH1F(("NumHits_layer_3_" + part_string ).c_str(), ("Number of hits (layer 3) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_4 = new TH1F(("NumHits_layer_4_" + part_string ).c_str(), ("Number of hits (layer 4) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_5 = new TH1F(("NumHits_layer_5_" + part_string ).c_str(), ("Number of hits (layer 5) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_6 = new TH1F(("NumHits_layer_6_" + part_string ).c_str(), ("Number of hits (layer 6) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_7 = new TH1F(("NumHits_layer_7_" + part_string ).c_str(), ("Number of hits (layer 7) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_8 = new TH1F(("NumHits_layer_8_" + part_string ).c_str(), ("Number of hits (layer 8) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_9 = new TH1F(("NumHits_layer_9_" + part_string ).c_str(), ("Number of hits (layer 9) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_10 = new TH1F(("NumHits_layer_10_" + part_string ).c_str(), ("Number of hits (layer 10) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_11 = new TH1F(("NumHits_layer_11_" + part_string ).c_str(), ("Number of hits (layer 11) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_12 = new TH1F(("NumHits_layer_12_" + part_string ).c_str(), ("Number of hits (layer 12) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_13 = new TH1F(("NumHits_layer_13_" + part_string ).c_str(), ("Number of hits (layer 13) " + part_string ).c_str(), 201, -0.5, 200.5);
    TH1F *h_nhit_layer_14 = new TH1F(("NumHits_layer_14_" + part_string ).c_str(), ("Number of hits (layer 14) " + part_string ).c_str(), 201, -0.5, 200.5);

    TH1F *h_nhit_layer_n_0 = new TH1F(("NumHits_layer_n_0_" + part_string ).c_str(), ("Number of hits normalized (layer 0) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_1 = new TH1F(("NumHits_layer_n_1_" + part_string ).c_str(), ("Number of hits normalized (layer 1) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_2 = new TH1F(("NumHits_layer_n_2_" + part_string ).c_str(), ("Number of hits normalized (layer 2) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_3 = new TH1F(("NumHits_layer_n_3_" + part_string ).c_str(), ("Number of hits normalized (layer 3) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_4 = new TH1F(("NumHits_layer_n_4_" + part_string ).c_str(), ("Number of hits normalized (layer 4) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_5 = new TH1F(("NumHits_layer_n_5_" + part_string ).c_str(), ("Number of hits normalized (layer 5) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_6 = new TH1F(("NumHits_layer_n_6_" + part_string ).c_str(), ("Number of hits normalized (layer 6) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_7 = new TH1F(("NumHits_layer_n_7_" + part_string ).c_str(), ("Number of hits normalized (layer 7) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_8 = new TH1F(("NumHits_layer_n_8_" + part_string ).c_str(), ("Number of hits normalized (layer 8) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_9 = new TH1F(("NumHits_layer_n_9_" + part_string ).c_str(), ("Number of hits normalized (layer 9) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_10 = new TH1F(("NumHits_layer_n_10_" + part_string ).c_str(), ("Number of hits normalized (layer 10) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_11 = new TH1F(("NumHits_layer_n_11_" + part_string ).c_str(), ("Number of hits normalized (layer 11) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_12 = new TH1F(("NumHits_layer_n_12_" + part_string ).c_str(), ("Number of hits normalized (layer 12) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_13 = new TH1F(("NumHits_layer_n_13_" + part_string ).c_str(), ("Number of hits normalized (layer 13) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_nhit_layer_n_14 = new TH1F(("NumHits_layer_n_14_" + part_string ).c_str(), ("Number of hits normalized (layer 14) " + part_string ).c_str(), 100, 0, 1);

    TH1F *h_weighte_layer_0 = new TH1F(("Weighte_layer_0_" + part_string ).c_str(), ("Weighted energy (layer 0) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_1 = new TH1F(("Weighte_layer_1_" + part_string ).c_str(), ("Weighted energy (layer 1) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_2 = new TH1F(("Weighte_layer_2_" + part_string ).c_str(), ("Weighted energy (layer 2) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_3 = new TH1F(("Weighte_layer_3_" + part_string ).c_str(), ("Weighted energy (layer 3) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_4 = new TH1F(("Weighte_layer_4_" + part_string ).c_str(), ("Weighted energy (layer 4) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_5 = new TH1F(("Weighte_layer_5_" + part_string ).c_str(), ("Weighted energy (layer 5) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_6 = new TH1F(("Weighte_layer_6_" + part_string ).c_str(), ("Weighted energy (layer 6) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_7 = new TH1F(("Weighte_layer_7_" + part_string ).c_str(), ("Weighted energy (layer 7) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_8 = new TH1F(("Weighte_layer_8_" + part_string ).c_str(), ("Weighted energy (layer 8) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_9 = new TH1F(("Weighte_layer_9_" + part_string ).c_str(), ("Weighted energy (layer 9) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_10 = new TH1F(("Weighte_layer_10_" + part_string ).c_str(), ("Weighted energy (layer 10) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_11 = new TH1F(("Weighte_layer_11_" + part_string ).c_str(), ("Weighted energy (layer 11) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_12 = new TH1F(("Weighte_layer_12_" + part_string ).c_str(), ("Weighted energy (layer 12) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_13 = new TH1F(("Weighte_layer_13_" + part_string ).c_str(), ("Weighted energy (layer 13) " + part_string ).c_str(), 1001, -0.5, 10000.5);
    TH1F *h_weighte_layer_14 = new TH1F(("Weighte_layer_14_" + part_string ).c_str(), ("Weighted energy (layer 14) " + part_string ).c_str(), 1001, -0.5, 10000.5);

    TH1F *h_weighte_layer_n_0 = new TH1F(("Weighte_layer_n_0_" + part_string ).c_str(), ("Weighted energy normalized (layer 0) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_1 = new TH1F(("Weighte_layer_n_1_" + part_string ).c_str(), ("Weighted energy normalized (layer 1) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_2 = new TH1F(("Weighte_layer_n_2_" + part_string ).c_str(), ("Weighted energy normalized (layer 2) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_3 = new TH1F(("Weighte_layer_n_3_" + part_string ).c_str(), ("Weighted energy normalized (layer 3) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_4 = new TH1F(("Weighte_layer_n_4_" + part_string ).c_str(), ("Weighted energy normalized (layer 4) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_5 = new TH1F(("Weighte_layer_n_5_" + part_string ).c_str(), ("Weighted energy normalized (layer 5) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_6 = new TH1F(("Weighte_layer_n_6_" + part_string ).c_str(), ("Weighted energy normalized (layer 6) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_7 = new TH1F(("Weighte_layer_n_7_" + part_string ).c_str(), ("Weighted energy normalized (layer 7) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_8 = new TH1F(("Weighte_layer_n_8_" + part_string ).c_str(), ("Weighted energy normalized (layer 8) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_9 = new TH1F(("Weighte_layer_n_9_" + part_string ).c_str(), ("Weighted energy normalized (layer 9) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_10 = new TH1F(("Weighte_layer_n_10_" + part_string ).c_str(), ("Weighted energy normalized (layer 10) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_11 = new TH1F(("Weighte_layer_n_11_" + part_string ).c_str(), ("Weighted energy normalized (layer 11) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_12 = new TH1F(("Weighte_layer_n_12_" + part_string ).c_str(), ("Weighted energy normalized (layer 12) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_13 = new TH1F(("Weighte_layer_n_13_" + part_string ).c_str(), ("Weighted energy normalized (layer 13) " + part_string ).c_str(), 100, 0, 1);
    TH1F *h_weighte_layer_n_14 = new TH1F(("Weighte_layer_n_14_" + part_string ).c_str(), ("Weighted energy normalized (layer 14) " + part_string ).c_str(), 100, 0, 1);



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
        	
	for (int i_event = 0; i_event < (int)nentries; i_event++) {
	  // Resolution
	  int nhit = 0;     
	  float sume = 0;   
	  float weighte = 0; 
	  float mol_value = 0.;
	  float bar_xyzr[4]; // 4 for (x,y,z,r)
	  
	  float nhit_layer_array[N_ECAL_LAYERS];
	  float nhit_layer_n_array[N_ECAL_LAYERS];

	  float weighte_layer_array[N_ECAL_LAYERS];
          float weighte_layer_n_array[N_ECAL_LAYERS];

	  float bar_layer_array[N_ECAL_LAYERS][3]; // 3 for (x,y,r)

	  tree->GetEntry(i_event);
	  if(i_event % 1000 == 0) cout << "Event " << to_string(i_event) << endl;
	  
	  get_res(nhit,	sume, weighte, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked);

	  h_nhit->Fill(nhit);        
	  h_sume->Fill(sume);         
	  h_weighte->Fill(weighte);     
	  
	  // Fill barycenter
	  barycenter(hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, bar_xyzr, hit_isMasked, masked);
	  
	  h_bar_x->Fill(bar_xyzr[0]);
	  h_bar_y->Fill(bar_xyzr[1]);
	  h_bar_z->Fill(bar_xyzr[2]);
	  h_bar_r->Fill(bar_xyzr[3]);
	  
	  // Fill Moliere radii histograms
	  mol_value = moliere(hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, hit_isMasked, masked);
	  h_mol->Fill(mol_value);
	  
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
          hits_layer(weighte_layer_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, false, "weight");
	  hits_layer(weighte_layer_n_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, true, "weight");
          
          h_weighte_layer_0->Fill(weighte_layer_array[0]);
          h_weighte_layer_1->Fill(weighte_layer_array[1]);
          h_weighte_layer_2->Fill(weighte_layer_array[2]);
          h_weighte_layer_3->Fill(weighte_layer_array[3]);
          h_weighte_layer_4->Fill(weighte_layer_array[4]);
          h_weighte_layer_5->Fill(weighte_layer_array[5]);
          h_weighte_layer_6->Fill(weighte_layer_array[6]);
          h_weighte_layer_7->Fill(weighte_layer_array[7]);
          h_weighte_layer_8->Fill(weighte_layer_array[8]);
          h_weighte_layer_9->Fill(weighte_layer_array[9]);
          h_weighte_layer_10->Fill(weighte_layer_array[10]);
          h_weighte_layer_11->Fill(weighte_layer_array[11]);
          h_weighte_layer_12->Fill(weighte_layer_array[12]);
          h_weighte_layer_13->Fill(weighte_layer_array[13]);
          h_weighte_layer_14->Fill(weighte_layer_array[14]);

	  h_weighte_layer_n_0->Fill(weighte_layer_n_array[0]);
          h_weighte_layer_n_1->Fill(weighte_layer_n_array[1]);
          h_weighte_layer_n_2->Fill(weighte_layer_n_array[2]);
          h_weighte_layer_n_3->Fill(weighte_layer_n_array[3]);
          h_weighte_layer_n_4->Fill(weighte_layer_n_array[4]);
          h_weighte_layer_n_5->Fill(weighte_layer_n_array[5]);
          h_weighte_layer_n_6->Fill(weighte_layer_n_array[6]);
          h_weighte_layer_n_7->Fill(weighte_layer_n_array[7]);
          h_weighte_layer_n_8->Fill(weighte_layer_n_array[8]);
          h_weighte_layer_n_9->Fill(weighte_layer_n_array[9]);
          h_weighte_layer_n_10->Fill(weighte_layer_n_array[10]);
          h_weighte_layer_n_11->Fill(weighte_layer_n_array[11]);
          h_weighte_layer_n_12->Fill(weighte_layer_n_array[12]);
          h_weighte_layer_n_13->Fill(weighte_layer_n_array[13]);
          h_weighte_layer_n_14->Fill(weighte_layer_n_array[14]);

	  // shower weight general parameters start
          float weighte_shower_maxvalue = 0.;
	  float weighte_shower_maxvalue_n = 0.;
          int weighte_ilayermax = -1;
          int weighte_ilayerstart = -1;
          int weighte_ilayerend = -1;
          int weighte_ilayerstart_10 = -1;
          int weighte_ilayerend_10 = -1;

	  float weighte_shower_averagevalue = weighte/N_ECAL_LAYERS;

	  shower_variables(weighte, weighte_layer_array, weighte_layer_n_array, weighte_shower_maxvalue, weighte_shower_maxvalue_n, weighte_ilayermax,
                           weighte_ilayerstart, weighte_ilayerstart_10, weighte_ilayerend, weighte_ilayerend_10, "weight");
	  //cout<<"weight shower variables: maxvalue, maxvalue_n,ilayermax, ilayerstart, ilayerstart(10%), ilayerend, ilayerend(10%)"<<endl;
	  //cout<<"weight shower variables: "<<weight<<" "<<weighte_shower_maxvalue<<" "<<weighte_shower_maxvalue_n<<" "<<weighte_ilayermax<<" "<<
          //  weighte_ilayerstart<<" "<<weighte_ilayerstart_10<<" "<<weighte_ilayerend<<" "<<weighte_ilayerend_10<<endl;

	  //shower weight general parameters finish	  
          h_shower_weighte_max_layer->Fill(weighte_ilayermax);
          h_shower_weighte_start_layer->Fill(weighte_ilayerstart);
          h_shower_weighte_end_layer->Fill(weighte_ilayerend);
          h_shower_weighte_start_10_layer->Fill(weighte_ilayerstart_10);
          h_shower_weighte_end_10_layer->Fill(weighte_ilayerend_10);
          h_shower_weighte_average->Fill(weighte_shower_averagevalue);
          h_shower_weighte_max->Fill(weighte_shower_maxvalue);
	  
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

	  // Filling the tree
	  b_nhit = nhit;
	  b_sume = sume;
	  b_weighte = weighte;
	  b_bar_x = bar_xyzr[0];
	  b_bar_y = bar_xyzr[1];
	  b_bar_z = bar_xyzr[2];
	  b_bar_r = bar_xyzr[3];
	  b_mol = mol_value;

	  b_nhit_layer_0 = nhit_layer_array[0];
	  b_nhit_layer_1 = nhit_layer_array[1];
	  b_nhit_layer_2 = nhit_layer_array[2];
	  b_nhit_layer_3 = nhit_layer_array[3];
	  b_nhit_layer_4 = nhit_layer_array[4];
	  b_nhit_layer_5 = nhit_layer_array[5];
	  b_nhit_layer_6 = nhit_layer_array[6];
	  b_nhit_layer_7 = nhit_layer_array[7];
	  b_nhit_layer_8 = nhit_layer_array[8];
	  b_nhit_layer_9 = nhit_layer_array[9];
	  b_nhit_layer_10 = nhit_layer_array[10];
	  b_nhit_layer_11 = nhit_layer_array[11];
	  b_nhit_layer_12 = nhit_layer_array[12];
	  b_nhit_layer_13 = nhit_layer_array[13];
	  b_nhit_layer_14 = nhit_layer_array[14];

	  b_nhit_layer_n_0 = nhit_layer_n_array[0];
	  b_nhit_layer_n_1 = nhit_layer_n_array[1];
	  b_nhit_layer_n_2 = nhit_layer_n_array[2];
	  b_nhit_layer_n_3 = nhit_layer_n_array[3];
	  b_nhit_layer_n_4 = nhit_layer_n_array[4];
	  b_nhit_layer_n_5 = nhit_layer_n_array[5];
	  b_nhit_layer_n_6 = nhit_layer_n_array[6];
	  b_nhit_layer_n_7 = nhit_layer_n_array[7];
	  b_nhit_layer_n_8 = nhit_layer_n_array[8];
	  b_nhit_layer_n_9 = nhit_layer_n_array[9];
	  b_nhit_layer_n_10 = nhit_layer_n_array[10];
	  b_nhit_layer_n_11 = nhit_layer_n_array[11];
	  b_nhit_layer_n_12 = nhit_layer_n_array[12];
	  b_nhit_layer_n_13 = nhit_layer_n_array[13];
	  b_nhit_layer_n_14 = nhit_layer_n_array[14];

	  b_shower_nhit_max_layer = nhit_ilayermax;
	  b_shower_nhit_start_layer = nhit_ilayerstart;
	  b_shower_nhit_end_layer = nhit_ilayerend;
	  b_shower_nhit_start_10_layer = nhit_ilayerstart_10;
	  b_shower_nhit_end_10_layer = nhit_ilayerend_10;
	  b_shower_nhit_average = nhit_shower_averagevalue;
	  b_shower_nhit_max = nhit_shower_maxvalue;

	  b_weighte_layer_0 = weighte_layer_array[0];
	  b_weighte_layer_1 = weighte_layer_array[1];
	  b_weighte_layer_2 = weighte_layer_array[2];
	  b_weighte_layer_3 = weighte_layer_array[3];
	  b_weighte_layer_4 = weighte_layer_array[4];
	  b_weighte_layer_5 = weighte_layer_array[5];
	  b_weighte_layer_6 = weighte_layer_array[6];
	  b_weighte_layer_7 = weighte_layer_array[7];
	  b_weighte_layer_8 = weighte_layer_array[8];
	  b_weighte_layer_9 = weighte_layer_array[9];
	  b_weighte_layer_10 = weighte_layer_array[10];
	  b_weighte_layer_11 = weighte_layer_array[11];
	  b_weighte_layer_12 = weighte_layer_array[12];
	  b_weighte_layer_13 = weighte_layer_array[13];
	  b_weighte_layer_14 = weighte_layer_array[14];

	  b_weighte_layer_n_0 = weighte_layer_n_array[0];
	  b_weighte_layer_n_1 = weighte_layer_n_array[1];
	  b_weighte_layer_n_2 = weighte_layer_n_array[2];
	  b_weighte_layer_n_3 = weighte_layer_n_array[3];
	  b_weighte_layer_n_4 = weighte_layer_n_array[4];
	  b_weighte_layer_n_5 = weighte_layer_n_array[5];
	  b_weighte_layer_n_6 = weighte_layer_n_array[6];
	  b_weighte_layer_n_7 = weighte_layer_n_array[7];
	  b_weighte_layer_n_8 = weighte_layer_n_array[8];
	  b_weighte_layer_n_9 = weighte_layer_n_array[9];
	  b_weighte_layer_n_10 = weighte_layer_n_array[10];
	  b_weighte_layer_n_11 = weighte_layer_n_array[11];
	  b_weighte_layer_n_12 = weighte_layer_n_array[12];
	  b_weighte_layer_n_13 = weighte_layer_n_array[13];
	  b_weighte_layer_n_14 = weighte_layer_n_array[14];

	  b_shower_weighte_max_layer = weighte_ilayermax;
	  b_shower_weighte_start_layer = weighte_ilayerstart;
	  b_shower_weighte_end_layer = weighte_ilayerend;
	  b_shower_weighte_start_10_layer = weighte_ilayerstart_10;
	  b_shower_weighte_end_10_layer = weighte_ilayerend_10;
	  b_shower_weighte_average = weighte_shower_averagevalue;
	  b_shower_weighte_max = weighte_shower_maxvalue;

	  b_bar_x_layer_0 = bar_layer_array[0][0];
	  b_bar_x_layer_1 = bar_layer_array[1][0];
	  b_bar_x_layer_2 = bar_layer_array[2][0];
	  b_bar_x_layer_3 = bar_layer_array[3][0];
	  b_bar_x_layer_4 = bar_layer_array[4][0];
	  b_bar_x_layer_5 = bar_layer_array[5][0];
	  b_bar_x_layer_6 = bar_layer_array[6][0];
	  b_bar_x_layer_7 = bar_layer_array[7][0];
	  b_bar_x_layer_8 = bar_layer_array[8][0];
	  b_bar_x_layer_9 = bar_layer_array[9][0];
	  b_bar_x_layer_10 = bar_layer_array[10][0];
	  b_bar_x_layer_11 = bar_layer_array[11][0];
	  b_bar_x_layer_12 = bar_layer_array[12][0];
	  b_bar_x_layer_13 = bar_layer_array[13][0];
	  b_bar_x_layer_14 = bar_layer_array[14][0];

	  b_bar_y_layer_0 = bar_layer_array[0][1];
	  b_bar_y_layer_1 = bar_layer_array[1][1];
	  b_bar_y_layer_2 = bar_layer_array[2][1];
	  b_bar_y_layer_3 = bar_layer_array[3][1];
	  b_bar_y_layer_4 = bar_layer_array[4][1];
	  b_bar_y_layer_5 = bar_layer_array[5][1];
	  b_bar_y_layer_6 = bar_layer_array[6][1];
	  b_bar_y_layer_7 = bar_layer_array[7][1];
	  b_bar_y_layer_8 = bar_layer_array[8][1];
	  b_bar_y_layer_9 = bar_layer_array[9][1];
	  b_bar_y_layer_10 = bar_layer_array[10][1];
	  b_bar_y_layer_11 = bar_layer_array[11][1];
	  b_bar_y_layer_12 = bar_layer_array[12][1];
	  b_bar_y_layer_13 = bar_layer_array[13][1];
	  b_bar_y_layer_14 = bar_layer_array[14][1];

	  b_bar_r_layer_0 = bar_layer_array[0][2];
          b_bar_r_layer_1 = bar_layer_array[1][2];
          b_bar_r_layer_2 = bar_layer_array[2][2];
          b_bar_r_layer_3 = bar_layer_array[3][2];
          b_bar_r_layer_4 = bar_layer_array[4][2];
          b_bar_r_layer_5 = bar_layer_array[5][2];
          b_bar_r_layer_6 = bar_layer_array[6][2];
          b_bar_r_layer_7 = bar_layer_array[7][2];
          b_bar_r_layer_8 = bar_layer_array[8][2];
          b_bar_r_layer_9 = bar_layer_array[9][2];
          b_bar_r_layer_10 = bar_layer_array[10][2];
          b_bar_r_layer_11 = bar_layer_array[11][2];
          b_bar_r_layer_12 = bar_layer_array[12][2];
          b_bar_r_layer_13 = bar_layer_array[13][2];
          b_bar_r_layer_14 = bar_layer_array[14][2];

	  outtree->Fill();
	  
	  hit_isMasked->clear();
	  hit_energy->clear();
	  hit_slab->clear();
	  
        }
    }
    // End of writting histos	
    //Write objects
    f.WriteTObject(outtree);
    f.WriteTObject(h_nhit);    
    f.WriteTObject(h_sume);    
    f.WriteTObject(h_weighte);  
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
    
    f.WriteTObject(h_weighte_layer_0);
    f.WriteTObject(h_weighte_layer_1);
    f.WriteTObject(h_weighte_layer_2);
    f.WriteTObject(h_weighte_layer_3);
    f.WriteTObject(h_weighte_layer_4);
    f.WriteTObject(h_weighte_layer_5);
    f.WriteTObject(h_weighte_layer_6);
    f.WriteTObject(h_weighte_layer_7);
    f.WriteTObject(h_weighte_layer_8);
    f.WriteTObject(h_weighte_layer_9);
    f.WriteTObject(h_weighte_layer_10);
    f.WriteTObject(h_weighte_layer_11);
    f.WriteTObject(h_weighte_layer_12);
    f.WriteTObject(h_weighte_layer_13);
    f.WriteTObject(h_weighte_layer_14);
    
    f.WriteTObject(h_weighte_layer_n_0);
    f.WriteTObject(h_weighte_layer_n_1);
    f.WriteTObject(h_weighte_layer_n_2);
    f.WriteTObject(h_weighte_layer_n_3);
    f.WriteTObject(h_weighte_layer_n_4);
    f.WriteTObject(h_weighte_layer_n_5);
    f.WriteTObject(h_weighte_layer_n_6);
    f.WriteTObject(h_weighte_layer_n_7);
    f.WriteTObject(h_weighte_layer_n_8);
    f.WriteTObject(h_weighte_layer_n_9);
    f.WriteTObject(h_weighte_layer_n_10);
    f.WriteTObject(h_weighte_layer_n_11);
    f.WriteTObject(h_weighte_layer_n_12);
    f.WriteTObject(h_weighte_layer_n_13);
    f.WriteTObject(h_weighte_layer_n_14);
    
    f.WriteTObject(h_shower_weighte_max_layer);
    f.WriteTObject(h_shower_weighte_start_layer);
    f.WriteTObject(h_shower_weighte_end_layer);
    f.WriteTObject(h_shower_weighte_start_10_layer);
    f.WriteTObject(h_shower_weighte_end_10_layer);
    f.WriteTObject(h_shower_weighte_average);
    f.WriteTObject(h_shower_weighte_max);
    
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
    
    f.Close();
	
}
