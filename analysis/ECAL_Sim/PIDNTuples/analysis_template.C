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

#define N_ENERGIES XENERGIESX //18
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
            weight += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
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
      sume_w_total += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
      //cout<<"|DEBUG(hit_layer)| sume total: "<<sume_total<<", sume_w_total: "<<sume_w_total<<endl;
      //cout<<"layer and thickness"<<hit_slab->at(j)<<" "<<W_thicknesses[hit_slab->at(j)]<<endl;
    }
  }
  
  for (int ilayer=0; ilayer < N_ECAL_LAYERS; ilayer++) {
    if(hit_slab->size() > 0){
      for( int j = 0; j < hit_energy->size(); j++ ) {
	if( masked && hit_isMasked->at(j) == 1 ) continue;
	if( hit_slab->at(j) == ilayer ) {
	  hit_count[ilayer] += 1;
	  sume[ilayer] += hit_energy->at(j);
	  sume_w[ilayer] += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
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

bool is_Shower(float entries, float array[N_ECAL_LAYERS]) {
  float threshold = 3.;
  float shower_maxvalue = 0.;
  bool isShower = false;


  if(entries > 0){
    for(int ilayer=0; ilayer<N_ECAL_LAYERS; ilayer++){
      float thislayer = array[ilayer];
      if(thislayer > shower_maxvalue){
        shower_maxvalue = thislayer;
      }
    }
    for(int ilayer=0; ilayer<N_ECAL_LAYERS-2; ilayer++){
      float thislayer = array[ilayer];
      float nextlayer = array[ilayer+1];
      float nextnextlayer = array[ilayer+2];
      if((thislayer > threshold) and (nextlayer > thislayer) and (nextnextlayer > nextlayer) and (shower_maxvalue > 5.)){
        isShower = true;
        break;
      }
    }
  }
  return isShower;

}

void shower_variables(float entries, float array[N_ECAL_LAYERS], float array_n[N_ECAL_LAYERS], float &shower_maxvalue, float &shower_maxvalue_n, int &ilayermax,
                      int &ilayerstart, int &ilayerstart_10, int &ilayerend, int &ilayerend_10, string count_type = "nhit", bool is_shower=false) {

  float percentage = 0.1;
  float threshold = 3.;
 
  if((entries > 0) and (is_shower == true)){
    for(int ilayer=0; ilayer<N_ECAL_LAYERS; ilayer++){
      float thislayer = array[ilayer];
      float thislayer_n = array_n[ilayer];
      if(thislayer > shower_maxvalue){
        shower_maxvalue = thislayer;
        shower_maxvalue_n = thislayer_n;
        ilayermax = ilayer;
      }
    }
    for(int ilayer=N_ECAL_LAYERS-1; ilayer>ilayermax; ilayer--){
      float thislayer = array[ilayer];
      if(thislayer > threshold){
        ilayerend = ilayer;
        break;
      }
    }
    for(int ilayer=N_ECAL_LAYERS-1; ilayer>ilayermax; ilayer--){
      float thislayer = array[ilayer];
      if((thislayer > threshold) and (thislayer > percentage*shower_maxvalue)) {
        ilayerend_10 = ilayer;
        break;
      }
    }
 
    for(int ilayer=0; ilayer<ilayermax; ilayer++){
      float thislayer = array[ilayer];
       if((thislayer > threshold)) {
        ilayerstart = ilayer;
        break;
      }
    }
    for(int ilayer=0; ilayer<ilayermax; ilayer++){
      float thislayer = array[ilayer];
       if((thislayer > threshold) && (thislayer > percentage*shower_maxvalue)) {
        ilayerstart_10 = ilayer;
        break;
      }
    }

  }
  return 0;
}

float MIP_Likeness(float nhits_layer[N_ECAL_LAYERS]) {
  float score = 0.;
  for(int i=0; i<N_ECAL_LAYERS; i++) {
    if(nhits_layer[i] == 0) score -= 1.;
    if(nhits_layer[i] > 0) score += 1./(nhits_layer[i]);
  }
  score = (score/N_ECAL_LAYERS + 1.)/2.;
  return score;
}

void is_interaction(float &ecal_int, int nhit_e) {
  ecal_int = 0.;
  if(nhit_e > 0) ecal_int = 1.;
  return 0;
}


float hits_max_distance(vector<int> *hit_slab, vector<float> * hit_x, vector<float> * hit_y,
			vector<int> * hit_isMasked, bool masked=false) {
  float max_distance = 0.;
  if(hit_slab->size() > 1){
    for (int i = 0; i < hit_slab->size(); i++){
      if (masked && hit_isMasked->at(i) == 1) continue;
      int layer1 = hit_slab->at(i);
      float hit1_x = hit_x->at(i);
      float hit1_y = hit_y->at(i);
      for (int j = i+1; j < hit_slab->size()-1; j++){
	if (masked && hit_isMasked->at(j) == 1) continue;
        int layer2 = hit_slab->at(j);
	if(layer2 != layer1) continue;
        float hit2_x = hit_x->at(j);
        float hit2_y = hit_y->at(j);
        float distance = pow(pow(hit2_x-hit1_x,2)+pow(hit2_y-hit1_y,2),0.5);
        if(distance>max_distance) max_distance = distance;
      }
    }
  }
  return max_distance;
}

struct hitpair
{
  float hit_rs;
  float hit_es;
};

bool CompareHitsR(const hitpair hit1, const hitpair hit2) { return hit1.hit_rs < hit2.hit_rs; }

float moliere(vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses,
		   vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z,
		   vector<int> * hit_isMasked, bool masked=false, float containment = 0.90, bool is_shower=false) {
  float mol_rad = 0.;
  float weighte = 0.;
  float wx = 0.; float wy = 0.; float wz = 0.;
  float r = 0.;

  if((hit_energy->size() > 0) and (is_shower==true)){
    vector<hitpair> hits_vector;

    for (int j = 0; j < hit_energy->size(); j++){
      if (masked && hit_isMasked->at(j) == 1) continue;
      weighte += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
      wx += hit_x->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
      wy += hit_y->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
    }

    float bary_x = wx / weighte;
    float bary_y = wy / weighte;
    if( (weighte == 0) or (weighte < 0) ){
      bary_x = 0.;
      bary_y = 0.;
    }

    for (int j = 0; j < hit_energy->size(); j++){
      if (masked && hit_isMasked->at(j) == 1) continue;
      r = pow(pow((hit_x->at(j) - bary_x) , 2) + pow((hit_y->at(j) - bary_y), 2), 0.5);
      hits_vector.push_back({static_cast<float>(r),static_cast<float>(hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5))});
    }

    float mol_e = 0.;
    int mol_i=0;
    if(hits_vector.size() > 0){
      std::sort(hits_vector.begin(), hits_vector.end(), CompareHitsR);
      for (int j = 0; j < hits_vector.size(); j++) {
        mol_e += hits_vector.at(j).hit_es;
        if (mol_e >= containment * weighte){
          mol_i=j;
          break;
        }
      }

      if(mol_i<0) mol_i=0;
      mol_rad=hits_vector.at(mol_i).hit_rs;
    }
  }
  return mol_rad;
}

void radius_layer(float mol_per_layer[N_ECAL_LAYERS], vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses,
                  vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z,
                  vector<int> * hit_isMasked, bool masked=false, float containment = 0.90, bool is_shower=false) {

  if((hit_energy->size() > 0) and (is_shower==true)){

    for(int ilayer = 0; ilayer < N_ECAL_LAYERS; ilayer++){
      float mol_rad = 0.;
      float weighte = 0.;
      float wx = 0.; float wy = 0.;
      float r = 0.;
      int nhit_layer = 0;
      struct hitcontent
      {
        float hit_rs;
        float hit_es;
      };
      vector<hitcontent> hits_vector;

      for (int j = 0; j < hit_energy->size(); j++){
        if (masked && hit_isMasked->at(j) == 1) continue;
        if(hit_slab->at(j) == ilayer) {
          weighte += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
          wx += hit_x->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
          wy += hit_y->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
          nhit_layer += 1;
        }
      }

      float bary_x = wx / weighte;
      float bary_y = wy / weighte;
      if( (weighte == 0) or (weighte < 0) ){
        bary_x = 0.;
        bary_y = 0.;
      }

      for (int j = 0; j < hit_energy->size(); j++){
        if (masked && hit_isMasked->at(j) == 1) continue;
        if (hit_slab->at(j) == ilayer) {
          r = pow(pow((hit_x->at(j) - bary_x) , 2) + pow((hit_y->at(j) - bary_y), 2), 0.5);
          hits_vector.push_back({static_cast<float>(r),static_cast<float>(hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5))});
        }
      }
      float mol_e = 0.;
      int mol_i=0;
      if(hits_vector.size() > 0){

	std::sort(hits_vector.begin(), hits_vector.end(),
                  [](const auto& i, const auto& j) { return i.hit_rs < j.hit_rs; } );

        for (int j = 0; j < hits_vector.size(); j++) {
          mol_e += hits_vector.at(j).hit_es;
          if (mol_e >= containment * weighte){
            mol_i=j;
            break;
          }
        }

        if(mol_i<0) mol_i=0.;
        mol_per_layer[ilayer] = hits_vector.at(mol_i).hit_rs;
      }
      if(nhit_layer < 3) mol_per_layer[ilayer] = 0.;
      if(mol_per_layer[ilayer] < 1) mol_per_layer[ilayer] = 0.;
      hits_vector.clear();
    }
  }
  /*
  for(int iout = 0; iout < N_ECAL_LAYERS; iout++) cout<<mol_per_layer[iout]<<" ";
  cout<<endl;
  */
  return 0;
}

void barycenter(vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses,
                vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z, float bar_xyzr[4],
                vector<int> * hit_isMasked, bool masked=false, bool is_shower=false) {

  float sume = 0.;
  float wx = 0.; float wy = 0.; float wz = 0.;
  float bary_x = 0., bary_y = 0., bary_z = 0.;
  // Removing shower condition
  //if((hit_energy->size() > 0) and (is_shower == true)){
  if(hit_energy->size() > 0){
    for (int j = 0; j < hit_energy->size(); j++){
      if (masked && hit_isMasked->at(j) == 1) continue;
      //cout<<"DBG "<<j<<endl;                                                                                                                                                                            
      sume += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
      wx += hit_x->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
      wy += hit_y->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
      wz += hit_z->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
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
                vector<float> * hit_z, vector<int> * hit_isMasked, bool masked=false, bool is_shower=false) {

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
  // Removing the shower condition
  // if((hit_energy->size() > 0) and (is_shower == true)){
  if(hit_energy->size() > 0){
    for (int ilayer=0; ilayer < N_ECAL_LAYERS; ilayer++) {
      if(hit_slab->size() > 0){
        for (int j = 0; j < hit_energy->size(); j++) {
          if (masked && hit_isMasked->at(j) == 1) continue;
          if (hit_slab->at(j) == 0) {
            sume_w[ilayer] += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
            wx[ilayer] += hit_x->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
            wy[ilayer] += hit_y->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
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
    double test_e[N_ENERGIES]={XENERGYSTRINGX};
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
    string base_path = "XBASEPATHX";
    // for (int j = 0; j < N_ENERGIES; j++) filenames[j] = base_path + "/CONF11/build/ECAL_QGSP_BERT_conf6_e-_" + to_string((int)round(energies[j])) +  "GeV_5kevt_build_masked.root";
    for (int j = 0; j < N_ENERGIES; j++){
      std::string str = std::to_string (energies[j]);
      str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
      str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
      cout<<str<<endl;
      filenames[j] = base_path + "output_LCIO2Build_TB2022-06_"+particle+"_" + str +  "GeV.root";
    }
    TString result_name = "resolution_"+particle+"_result.root" ;
    
    TFile f(result_name, "recreate");
    TTree *outtree = new TTree("ntp","NTuples");

    // branches definitions
    // We will save both the branches and the histos
    // Values
    Float_t b_ecal_interaction;
    Float_t b_nhit, b_sume, b_weighte, b_bar_x, b_bar_y, b_bar_z, b_bar_r;
    Float_t b_mol;
    Float_t b_MIP_Likeness;
    Float_t b_hits_max_distance;
    Float_t b_radius90_layer_0, b_radius90_layer_1, b_radius90_layer_2, b_radius90_layer_3, b_radius90_layer_4, b_radius90_layer_5, b_radius90_layer_6, b_radius90_layer_7, b_radius90_layer_8, b_radius90_layer_9, b_radius90_layer_10, b_radius90_layer_11, b_radius90_layer_12, b_radius90_layer_13, b_radius90_layer_14;
    Float_t b_bar_x_layer_0, b_bar_x_layer_1, b_bar_x_layer_2, b_bar_x_layer_3, b_bar_x_layer_4, b_bar_x_layer_5, b_bar_x_layer_6, b_bar_x_layer_7, b_bar_x_layer_8, b_bar_x_layer_9, b_bar_x_layer_10, b_bar_x_layer_11, b_bar_x_layer_12, b_bar_x_layer_13, b_bar_x_layer_14;
    Float_t b_bar_y_layer_0, b_bar_y_layer_1, b_bar_y_layer_2, b_bar_y_layer_3, b_bar_y_layer_4, b_bar_y_layer_5, b_bar_y_layer_6, b_bar_y_layer_7, b_bar_y_layer_8, b_bar_y_layer_9, b_bar_y_layer_10, b_bar_y_layer_11, b_bar_y_layer_12, b_bar_y_layer_13, b_bar_y_layer_14;
    Float_t b_bar_r_layer_0, b_bar_r_layer_1, b_bar_r_layer_2, b_bar_r_layer_3, b_bar_r_layer_4, b_bar_r_layer_5, b_bar_r_layer_6, b_bar_r_layer_7, b_bar_r_layer_8, b_bar_r_layer_9, b_bar_r_layer_10, b_bar_r_layer_11, b_bar_r_layer_12, b_bar_r_layer_13, b_bar_r_layer_14;
    Float_t b_shower_nhit_max_layer, b_shower_nhit_start_layer, b_shower_nhit_end_layer, b_shower_nhit_start_10_layer, b_shower_nhit_end_10_layer, b_shower_nhit_average, b_shower_nhit_max;
    Float_t b_shower_sume_max_layer, b_shower_sume_start_layer, b_shower_sume_end_layer, b_shower_sume_start_10_layer, b_shower_sume_end_10_layer, b_shower_sume_average, b_shower_sume_max;
    Float_t b_shower_weighte_max_layer, b_shower_weighte_start_layer, b_shower_weighte_end_layer, b_shower_weighte_start_10_layer, b_shower_weighte_end_10_layer, b_shower_weighte_average, b_shower_weighte_max;
    Float_t b_nhit_layer_0, b_nhit_layer_1, b_nhit_layer_2, b_nhit_layer_3, b_nhit_layer_4, b_nhit_layer_5, b_nhit_layer_6, b_nhit_layer_7, b_nhit_layer_8, b_nhit_layer_9, b_nhit_layer_10, b_nhit_layer_11, b_nhit_layer_12, b_nhit_layer_13, b_nhit_layer_14;
    Float_t b_nhit_layer_n_0, b_nhit_layer_n_1, b_nhit_layer_n_2, b_nhit_layer_n_3, b_nhit_layer_n_4, b_nhit_layer_n_5, b_nhit_layer_n_6, b_nhit_layer_n_7, b_nhit_layer_n_8, b_nhit_layer_n_9, b_nhit_layer_n_10, b_nhit_layer_n_11, b_nhit_layer_n_12, b_nhit_layer_n_13, b_nhit_layer_n_14;
    Float_t b_weighte_layer_0, b_weighte_layer_1, b_weighte_layer_2, b_weighte_layer_3, b_weighte_layer_4, b_weighte_layer_5, b_weighte_layer_6, b_weighte_layer_7, b_weighte_layer_8, b_weighte_layer_9, b_weighte_layer_10, b_weighte_layer_11, b_weighte_layer_12, b_weighte_layer_13, b_weighte_layer_14;
    Float_t b_weighte_layer_n_0, b_weighte_layer_n_1, b_weighte_layer_n_2, b_weighte_layer_n_3, b_weighte_layer_n_4, b_weighte_layer_n_5, b_weighte_layer_n_6, b_weighte_layer_n_7, b_weighte_layer_n_8, b_weighte_layer_n_9, b_weighte_layer_n_10, b_weighte_layer_n_11, b_weighte_layer_n_12, b_weighte_layer_n_13, b_weighte_layer_n_14;    
    Float_t b_sume_layer_0, b_sume_layer_1, b_sume_layer_2, b_sume_layer_3, b_sume_layer_4, b_sume_layer_5, b_sume_layer_6, b_sume_layer_7, b_sume_layer_8, b_sume_layer_9, b_sume_layer_10, b_sume_layer_11, b_sume_layer_12, b_sume_layer_13, b_sume_layer_14;
    Float_t b_sume_layer_n_0, b_sume_layer_n_1, b_sume_layer_n_2, b_sume_layer_n_3, b_sume_layer_n_4, b_sume_layer_n_5, b_sume_layer_n_6, b_sume_layer_n_7, b_sume_layer_n_8, b_sume_layer_n_9, b_sume_layer_n_10, b_sume_layer_n_11, b_sume_layer_n_12, b_sume_layer_n_13, b_sume_layer_n_14;

    // Adding the branches
    outtree->Branch("ecal_interaction",&b_ecal_interaction,"b_ecal_interaction/F");
    outtree->Branch("nhit",&b_nhit,"b_nhit/F");
    outtree->Branch("sume",&b_sume,"b_sume/F");
    outtree->Branch("weighte",&b_weighte,"b_weighte/F");
    outtree->Branch("mol",&b_mol,"mol/F");
    outtree->Branch("MIP_Likeness",&b_MIP_Likeness,"MIP_Likeness/F");
    outtree->Branch("hits_max_distance",&b_hits_max_distance,"b_hits_max_distance/F");
    outtree->Branch("radius90_layer_0",&b_radius90_layer_0,"b_radius90_layer_0/F");
    outtree->Branch("radius90_layer_1",&b_radius90_layer_1,"b_radius90_layer_1/F");
    outtree->Branch("radius90_layer_2",&b_radius90_layer_2,"b_radius90_layer_2/F");
    outtree->Branch("radius90_layer_3",&b_radius90_layer_3,"b_radius90_layer_3/F");
    outtree->Branch("radius90_layer_4",&b_radius90_layer_4,"b_radius90_layer_4/F");
    outtree->Branch("radius90_layer_5",&b_radius90_layer_5,"b_radius90_layer_5/F");
    outtree->Branch("radius90_layer_6",&b_radius90_layer_6,"b_radius90_layer_6/F");
    outtree->Branch("radius90_layer_7",&b_radius90_layer_7,"b_radius90_layer_7/F");
    outtree->Branch("radius90_layer_8",&b_radius90_layer_8,"b_radius90_layer_8/F");
    outtree->Branch("radius90_layer_9",&b_radius90_layer_9,"b_radius90_layer_9/F");
    outtree->Branch("radius90_layer_10",&b_radius90_layer_10,"b_radius90_layer_10/F");
    outtree->Branch("radius90_layer_11",&b_radius90_layer_11,"b_radius90_layer_11/F");
    outtree->Branch("radius90_layer_12",&b_radius90_layer_12,"b_radius90_layer_12/F");
    outtree->Branch("radius90_layer_13",&b_radius90_layer_13,"b_radius90_layer_13/F");
    outtree->Branch("radius90_layer_14",&b_radius90_layer_14,"b_radius90_layer_14/F");
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
    outtree->Branch("shower_sume_max_layer",&b_shower_sume_max_layer,"b_shower_sume_max_layer/F");
    outtree->Branch("shower_sume_start_layer",&b_shower_sume_start_layer,"b_shower_sume_start_layer/F");
    outtree->Branch("shower_sume_end_layer",&b_shower_sume_end_layer,"b_shower_sume_end_layer/F");
    outtree->Branch("shower_sume_start_10_layer",&b_shower_sume_start_10_layer,"b_shower_sume_start_10_layer/F");
    outtree->Branch("shower_sume_end_10_layer",&b_shower_sume_end_10_layer,"b_shower_sume_end_10_layer/F");
    outtree->Branch("shower_sume_average",&b_shower_sume_average,"b_shower_sume_average/F");
    outtree->Branch("shower_sume_max",&b_shower_sume_max,"b_shower_sume_max/F");
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
    outtree->Branch("sume_layer_0",&b_sume_layer_0,"b_sume_layer_0/F");
    outtree->Branch("sume_layer_1",&b_sume_layer_1,"b_sume_layer_1/F");
    outtree->Branch("sume_layer_2",&b_sume_layer_2,"b_sume_layer_2/F");
    outtree->Branch("sume_layer_3",&b_sume_layer_3,"b_sume_layer_3/F");
    outtree->Branch("sume_layer_4",&b_sume_layer_4,"b_sume_layer_4/F");
    outtree->Branch("sume_layer_5",&b_sume_layer_5,"b_sume_layer_5/F");
    outtree->Branch("sume_layer_6",&b_sume_layer_6,"b_sume_layer_6/F");
    outtree->Branch("sume_layer_7",&b_sume_layer_7,"b_sume_layer_7/F");
    outtree->Branch("sume_layer_8",&b_sume_layer_8,"b_sume_layer_8/F");
    outtree->Branch("sume_layer_9",&b_sume_layer_9,"b_sume_layer_9/F");
    outtree->Branch("sume_layer_10",&b_sume_layer_10,"b_sume_layer_10/F");
    outtree->Branch("sume_layer_11",&b_sume_layer_11,"b_sume_layer_11/F");
    outtree->Branch("sume_layer_12",&b_sume_layer_12,"b_sume_layer_12/F");
    outtree->Branch("sume_layer_13",&b_sume_layer_13,"b_sume_layer_13/F");
    outtree->Branch("sume_layer_14",&b_sume_layer_14,"b_sume_layer_14/F");
    outtree->Branch("sume_layer_n_0",&b_sume_layer_n_0,"b_sume_layer_n_0/F");
    outtree->Branch("sume_layer_n_1",&b_sume_layer_n_1,"b_sume_layer_n_1/F");
    outtree->Branch("sume_layer_n_2",&b_sume_layer_n_2,"b_sume_layer_n_2/F");
    outtree->Branch("sume_layer_n_3",&b_sume_layer_n_3,"b_sume_layer_n_3/F");
    outtree->Branch("sume_layer_n_4",&b_sume_layer_n_4,"b_sume_layer_n_4/F");
    outtree->Branch("sume_layer_n_5",&b_sume_layer_n_5,"b_sume_layer_n_5/F");
    outtree->Branch("sume_layer_n_6",&b_sume_layer_n_6,"b_sume_layer_n_6/F");
    outtree->Branch("sume_layer_n_7",&b_sume_layer_n_7,"b_sume_layer_n_7/F");
    outtree->Branch("sume_layer_n_8",&b_sume_layer_n_8,"b_sume_layer_n_8/F");
    outtree->Branch("sume_layer_n_9",&b_sume_layer_n_9,"b_sume_layer_n_9/F");
    outtree->Branch("sume_layer_n_10",&b_sume_layer_n_10,"b_sume_layer_n_10/F");
    outtree->Branch("sume_layer_n_11",&b_sume_layer_n_11,"b_sume_layer_n_11/F");
    outtree->Branch("sume_layer_n_12",&b_sume_layer_n_12,"b_sume_layer_n_12/F");
    outtree->Branch("sume_layer_n_13",&b_sume_layer_n_13,"b_sume_layer_n_13/F");
    outtree->Branch("sume_layer_n_14",&b_sume_layer_n_14,"b_sume_layer_n_14/F");
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
    outtree->Branch("sume_layer_0",&b_sume_layer_0,"b_sume_layer_0/F");
    outtree->Branch("sume_layer_1",&b_sume_layer_1,"b_sume_layer_1/F");
    outtree->Branch("sume_layer_2",&b_sume_layer_2,"b_sume_layer_2/F");
    outtree->Branch("sume_layer_3",&b_sume_layer_3,"b_sume_layer_3/F");
    outtree->Branch("sume_layer_4",&b_sume_layer_4,"b_sume_layer_4/F");
    outtree->Branch("sume_layer_5",&b_sume_layer_5,"b_sume_layer_5/F");
    outtree->Branch("sume_layer_6",&b_sume_layer_6,"b_sume_layer_6/F");
    outtree->Branch("sume_layer_7",&b_sume_layer_7,"b_sume_layer_7/F");
    outtree->Branch("sume_layer_8",&b_sume_layer_8,"b_sume_layer_8/F");
    outtree->Branch("sume_layer_9",&b_sume_layer_9,"b_sume_layer_9/F");
    outtree->Branch("sume_layer_10",&b_sume_layer_10,"b_sume_layer_10/F");
    outtree->Branch("sume_layer_11",&b_sume_layer_11,"b_sume_layer_11/F");
    outtree->Branch("sume_layer_12",&b_sume_layer_12,"b_sume_layer_12/F");
    outtree->Branch("sume_layer_13",&b_sume_layer_13,"b_sume_layer_13/F");
    outtree->Branch("sume_layer_14",&b_sume_layer_14,"b_sume_layer_14/F");
    outtree->Branch("sume_layer_n_0",&b_sume_layer_n_0,"b_sume_layer_n_0/F");
    outtree->Branch("sume_layer_n_1",&b_sume_layer_n_1,"b_sume_layer_n_1/F");
    outtree->Branch("sume_layer_n_2",&b_sume_layer_n_2,"b_sume_layer_n_2/F");
    outtree->Branch("sume_layer_n_3",&b_sume_layer_n_3,"b_sume_layer_n_3/F");
    outtree->Branch("sume_layer_n_4",&b_sume_layer_n_4,"b_sume_layer_n_4/F");
    outtree->Branch("sume_layer_n_5",&b_sume_layer_n_5,"b_sume_layer_n_5/F");
    outtree->Branch("sume_layer_n_6",&b_sume_layer_n_6,"b_sume_layer_n_6/F");
    outtree->Branch("sume_layer_n_7",&b_sume_layer_n_7,"b_sume_layer_n_7/F");
    outtree->Branch("sume_layer_n_8",&b_sume_layer_n_8,"b_sume_layer_n_8/F");
    outtree->Branch("sume_layer_n_9",&b_sume_layer_n_9,"b_sume_layer_n_9/F");
    outtree->Branch("sume_layer_n_10",&b_sume_layer_n_10,"b_sume_layer_n_10/F");
    outtree->Branch("sume_layer_n_11",&b_sume_layer_n_11,"b_sume_layer_n_11/F");
    outtree->Branch("sume_layer_n_12",&b_sume_layer_n_12,"b_sume_layer_n_12/F");
    outtree->Branch("sume_layer_n_13",&b_sume_layer_n_13,"b_sume_layer_n_13/F");
    outtree->Branch("sume_layer_n_14",&b_sume_layer_n_14,"b_sume_layer_n_14/F");

    // Histo definitions, before the energy loop
    
    // Resolution histos
    string part_string = particle;

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

	  float sume_layer_array[N_ECAL_LAYERS];
          float sume_layer_n_array[N_ECAL_LAYERS];

	  float weighte_layer_array[N_ECAL_LAYERS];
          float weighte_layer_n_array[N_ECAL_LAYERS];

	  float bar_layer_array[N_ECAL_LAYERS][3]; // 3 for (x,y,r)
	  float radius90_layer_array[N_ECAL_LAYERS];

	  for(int i = 0; i<N_ECAL_LAYERS;i++) {
            nhit_layer_array[i] = 0.;
            nhit_layer_n_array[i] = 0.;
            sume_layer_array[i] = 0.;
            sume_layer_n_array[i] = 0.;
            weighte_layer_array[i] = 0.;
            weighte_layer_n_array[i] = 0.;
            bar_layer_array[i][0] = 0.;
            bar_layer_array[i][1] = 0.;
            bar_layer_array[i][2] = 0.;
            radius90_layer_array[i] = 0.;
          }
	  	  	  
	  tree->GetEntry(i_event);
	  if(i_event % 1000 == 0) cout << "Event " << to_string(i_event) << endl;
	  
	  get_res(nhit,	sume, weighte, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked);
	  
	  // Fill barycenter
	  barycenter(hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, bar_xyzr, hit_isMasked, masked);
	  
	  // Fill shower profile 
	  // Nhit
	  hits_layer(nhit_layer_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, false, "nhit");
	  hits_layer(nhit_layer_n_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, true, "nhit");
	  
	  // Fill MIP-Likeness
	  float MIP_Likeness_value = MIP_Likeness(nhit_layer_array);

	  // Sume
          hits_layer(sume_layer_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, false, "sume");
          hits_layer(sume_layer_n_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, true, "sume");

	  // Weighted energy
          hits_layer(weighte_layer_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, false, "weight");
          hits_layer(weighte_layer_n_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, true, "weight");

	  bool shower_bool = is_Shower(sume, sume_layer_array);

	  is_interaction(b_ecal_interaction, nhit);
	  
	  
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
			   nhit_ilayerstart, nhit_ilayerstart_10, nhit_ilayerend, nhit_ilayerend_10, "nhit", shower_bool);
	  //cout<<"nhit shower variables: maxvalue, maxvalue_n,ilayermax, ilayerstart, ilayerstart(10%), ilayerend, ilayerend(10%)"<<endl;
	  //cout<<"nhit shower variables: "<<nhit<<" "<<nhit_shower_maxvalue<<" "<<nhit_shower_maxvalue_n<<" "<<nhit_ilayermax<<" "<<
	  //  nhit_ilayerstart<<" "<<nhit_ilayerstart_10<<" "<<nhit_ilayerend<<" "<<nhit_ilayerend_10<<endl;

	  //shower nhit general parameters finish
	  
          // shower sume general parameters start
          float sume_shower_maxvalue = 0.;
          float sume_shower_maxvalue_n = 0.;
          int sume_ilayermax = -1;
          int sume_ilayerstart = -1;
          int sume_ilayerend = -1;
          int sume_ilayerstart_10 = -1;
          int sume_ilayerend_10 = -1;

          float sume_shower_averagevalue = sume/N_ECAL_LAYERS;

          shower_variables(sume, sume_layer_array, sume_layer_n_array, sume_shower_maxvalue, sume_shower_maxvalue_n, sume_ilayermax,
                           sume_ilayerstart, sume_ilayerstart_10, sume_ilayerend, sume_ilayerend_10, "sume", shower_bool);
         
          //shower sume general parameters finish

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
                           weighte_ilayerstart, weighte_ilayerstart_10, weighte_ilayerend, weighte_ilayerend_10, "weight", shower_bool);
	  //cout<<"weight shower variables: maxvalue, maxvalue_n,ilayermax, ilayerstart, ilayerstart(10%), ilayerend, ilayerend(10%)"<<endl;
	  //cout<<"weight shower variables: "<<weight<<" "<<weighte_shower_maxvalue<<" "<<weighte_shower_maxvalue_n<<" "<<weighte_ilayermax<<" "<<
          //  weighte_ilayerstart<<" "<<weighte_ilayerstart_10<<" "<<weighte_ilayerend<<" "<<weighte_ilayerend_10<<endl;

	  //shower weight general parameters finish	  

	  // Fill Moliere radii histograms 
          mol_value = moliere(hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, hit_isMasked, masked, 0.9, shower_bool);
          radius_layer(radius90_layer_array, hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, hit_isMasked, masked, 0.9, shower_bool);
          	  
	  // Barycenter x and y
	  bary_layer(bar_layer_array, hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, hit_isMasked, masked);
	  
	  // Hits max distance in the same layer
	  float hits_max_distance_value = hits_max_distance(hit_slab, hit_x, hit_y, hit_isMasked, masked);
	  
	  // Filling the tree
	  b_nhit = nhit;
	  b_sume = sume;
	  b_weighte = weighte;
	  b_bar_x = bar_xyzr[0];
	  b_bar_y = bar_xyzr[1];
	  b_bar_z = bar_xyzr[2];
	  b_bar_r = bar_xyzr[3];
	  b_mol = mol_value;
	  b_MIP_Likeness = MIP_Likeness_value;
	  b_hits_max_distance = hits_max_distance_value;

	  b_radius90_layer_0 = radius90_layer_array[0];
          b_radius90_layer_1 = radius90_layer_array[1];
          b_radius90_layer_2 = radius90_layer_array[2];
          b_radius90_layer_3 = radius90_layer_array[3];
          b_radius90_layer_4 = radius90_layer_array[4];
          b_radius90_layer_5 = radius90_layer_array[5];
          b_radius90_layer_6 = radius90_layer_array[6];
          b_radius90_layer_7 = radius90_layer_array[7];
          b_radius90_layer_8 = radius90_layer_array[8];
          b_radius90_layer_9 = radius90_layer_array[9];
          b_radius90_layer_10 = radius90_layer_array[10];
          b_radius90_layer_11 = radius90_layer_array[11];
          b_radius90_layer_12 = radius90_layer_array[12];
          b_radius90_layer_13 = radius90_layer_array[13];
          b_radius90_layer_14 = radius90_layer_array[14];

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

	  b_sume_layer_0 = sume_layer_array[0];
          b_sume_layer_1 = sume_layer_array[1];
          b_sume_layer_2 = sume_layer_array[2];
          b_sume_layer_3 = sume_layer_array[3];
          b_sume_layer_4 = sume_layer_array[4];
          b_sume_layer_5 = sume_layer_array[5];
          b_sume_layer_6 = sume_layer_array[6];
          b_sume_layer_7 = sume_layer_array[7];
          b_sume_layer_8 = sume_layer_array[8];
          b_sume_layer_9 = sume_layer_array[9];
          b_sume_layer_10 = sume_layer_array[10];
          b_sume_layer_11 = sume_layer_array[11];
          b_sume_layer_12 = sume_layer_array[12];
          b_sume_layer_13 = sume_layer_array[13];
          b_sume_layer_14 = sume_layer_array[14];

          b_sume_layer_n_0 = sume_layer_n_array[0];
          b_sume_layer_n_1 = sume_layer_n_array[1];
          b_sume_layer_n_2 = sume_layer_n_array[2];
          b_sume_layer_n_3 = sume_layer_n_array[3];
          b_sume_layer_n_4 = sume_layer_n_array[4];
          b_sume_layer_n_5 = sume_layer_n_array[5];
          b_sume_layer_n_6 = sume_layer_n_array[6];
          b_sume_layer_n_7 = sume_layer_n_array[7];
          b_sume_layer_n_8 = sume_layer_n_array[8];
          b_sume_layer_n_9 = sume_layer_n_array[9];
          b_sume_layer_n_10 = sume_layer_n_array[10];
          b_sume_layer_n_11 = sume_layer_n_array[11];
          b_sume_layer_n_12 = sume_layer_n_array[12];
          b_sume_layer_n_13 = sume_layer_n_array[13];
          b_sume_layer_n_14 = sume_layer_n_array[14];

          b_shower_sume_max_layer = sume_ilayermax;
          b_shower_sume_start_layer = sume_ilayerstart;
          b_shower_sume_end_layer = sume_ilayerend;
          b_shower_sume_start_10_layer = sume_ilayerstart_10;
          b_shower_sume_end_10_layer = sume_ilayerend_10;
          b_shower_sume_average = sume_shower_averagevalue;
          b_shower_sume_max = sume_shower_maxvalue;

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
    //Write objects
    f.WriteTObject(outtree);
        
    f.Close();
	
}
