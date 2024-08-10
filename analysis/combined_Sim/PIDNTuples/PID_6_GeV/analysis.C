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
#define N_HCAL_LAYERS 41 //pre+38+bad+tokyo

#define N_TOTAL_LAYERS 56

void get_energy_ecal(int &nhit, float &sume, float &weight, vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses, vector<int> *hit_isMasked, bool &masked) {
  if(hit_energy->size() > 0){    
    //cout<<W_thicknesses.Min()<<endl;
    // First option: use the minimum of the used
    // Second option: use the paper as reference 0.4X0, X0=3.5mm
    // It was weight_masked += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
    // New version: Divided by W X0
    
   for (int j = 0; j < hit_energy->size(); j++) {
     if( masked && hit_isMasked->at(j) == 1 ) continue;
            nhit += 1;
            sume += hit_energy->at(j);
            weight += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
   }
  }
  return;
}

void get_energy_hcal(int &nhit, float &sume, float &weight, vector<float> *hit_energy, vector<int> *hit_slab, TVectorD S_thicknesses, vector<int> *hit_isMasked, bool &masked) {
  // Divided by X0 of Steel: 17.618mm, read from materialScan
  // Steel thickness preshower: 0.5mm
  // Steel thickness per layer (next 38): (0.5 + 17.2 + 0.5) = 18.2mm
  // Steel thickness before bad layer: (0.5 + 17.2 + 17.2 + 0.5) = 35.4mm
  // Steel thickness before tokyo layer: (0.5 + 17.2 + 17.2 + 0.5) = 35.4mm
  
  if(hit_energy->size() > 0){
    for (int j = 0; j < hit_energy->size(); j++) {
      if( masked && hit_isMasked->at(j) == 1 ) continue;
      nhit += 1;
      sume += hit_energy->at(j);
      weight += hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
    }
  }
  return;
}


void hits_layer_ecal(float hlv[N_ECAL_LAYERS], vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses, vector<int> *hit_isMasked, bool masked=false, bool normalized=false, string count_type="nhit"){
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

void hits_layer_hcal(float hlv[N_HCAL_LAYERS], vector<float> * hit_energy, vector<int> *hit_slab, TVectorD S_thicknesses, vector<int> *hit_isMasked, bool masked=false, bool normalized=false, string count_type="nhit"){
  float hit_count[N_HCAL_LAYERS];
  float sume[N_HCAL_LAYERS];
  float sume_w[N_HCAL_LAYERS];

  float sume_total = 1.;
  float sume_w_total = 1.;

  float weight = 1.;

  for (int ilayer=0; ilayer < N_HCAL_LAYERS; ilayer++){
    hit_count[ilayer] = 0. ;
    sume[ilayer] = 0. ;
    sume_w[ilayer] = 0. ;
  }

  if(hit_energy->size() > 0){
    for (int j = 0; j < hit_energy->size(); j++) {
      if( masked && hit_isMasked->at(j) == 1 ) continue;
      sume_total += hit_energy->at(j);
      sume_w_total += hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
    }
  }

  for (int ilayer=0; ilayer < N_HCAL_LAYERS; ilayer++) {
    if(hit_slab->size() > 0){
      for( int j = 0; j < hit_energy->size(); j++ ) {
        if( masked && hit_isMasked->at(j) == 1 ) continue;
        if( hit_slab->at(j) == ilayer ) {
          hit_count[ilayer] += 1;
          sume[ilayer] += hit_energy->at(j);
	  sume_w[ilayer] += hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
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

return 0;
}


bool is_Shower_ecal(float entries, float array[N_ECAL_LAYERS]) {
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

bool is_Shower_hcal(float entries, float array[N_HCAL_LAYERS]) {
  float threshold = 3.;
  float shower_maxvalue = 0.;
  bool isShower = false;
  if(entries > 0){
    for(int ilayer=0; ilayer<N_HCAL_LAYERS; ilayer++){
      float thislayer = array[ilayer];
      if(thislayer > shower_maxvalue){
        shower_maxvalue = thislayer;
      }
    }
    for(int ilayer=0; ilayer<N_HCAL_LAYERS-2; ilayer++){
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

void shower_variables_ecal(float entries, float array[N_ECAL_LAYERS], float array_n[N_ECAL_LAYERS], float &shower_maxvalue, float &shower_maxvalue_n, int &ilayermax,
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

void shower_variables_hcal(float entries, float array[N_HCAL_LAYERS], float array_n[N_HCAL_LAYERS], float &shower_maxvalue, float &shower_maxvalue_n, int &ilayermax,
			   int &ilayerstart, int &ilayerstart_10, int &ilayerend, int &ilayerend_10, string count_type = "nhit", bool is_shower=false) {

  float percentage = 0.1;
  float threshold = 3.;

  if((entries > 0) and (is_shower == true)){
    for(int ilayer=0; ilayer<N_HCAL_LAYERS; ilayer++){
      float thislayer = array[ilayer];
      float thislayer_n = array_n[ilayer];
      if(thislayer > shower_maxvalue){
        shower_maxvalue = thislayer;
        shower_maxvalue_n = thislayer_n;
        ilayermax = ilayer;
      }
    }
    for(int ilayer=N_HCAL_LAYERS-1; ilayer>ilayermax; ilayer--){
      float thislayer = array[ilayer];
      if(thislayer > threshold){
        ilayerend = ilayer;
        break;
      }
    }
    for(int ilayer=N_HCAL_LAYERS-1; ilayer>ilayermax; ilayer--){
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

void shower_variables_total(float entries, float array[N_TOTAL_LAYERS], float array_n[N_TOTAL_LAYERS], float &shower_maxvalue, float &shower_maxvalue_n, int &ilayermax,
                           int &ilayerstart, int &ilayerstart_10, int &ilayerend, int &ilayerend_10, string count_type = "nhit", bool is_shower=false) {

  float percentage = 0.1;
  float threshold = 3.;

  if((entries > 0) and (is_shower == true)){
    for(int ilayer=0; ilayer<N_HCAL_LAYERS; ilayer++){
      float thislayer = array[ilayer];
      float thislayer_n = array_n[ilayer];
      if(thislayer > shower_maxvalue){
        shower_maxvalue = thislayer;
        shower_maxvalue_n = thislayer_n;
        ilayermax = ilayer;
      }
    }
    for(int ilayer=N_HCAL_LAYERS-1; ilayer>ilayermax; ilayer--){
      float thislayer = array[ilayer];
      if(thislayer > threshold){
        ilayerend = ilayer;
        break;
      }
    }
    for(int ilayer=N_HCAL_LAYERS-1; ilayer>ilayermax; ilayer--){
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

float MIP_Likeness_ecal(float nhits_layer[N_ECAL_LAYERS]) {
  float score = 0.;
  for(int i=0; i<N_ECAL_LAYERS; i++) {
    if(nhits_layer[i] == 0) score -= 1.;
    if(nhits_layer[i] > 0) score += 1./(nhits_layer[i]);
  }
  score = (score/N_ECAL_LAYERS + 1.)/2.;
  return score;
}

float MIP_Likeness_hcal(float nhits_layer[N_HCAL_LAYERS]) {
  float score = 0.;
  for(int i=0; i<N_HCAL_LAYERS; i++) {
    if(nhits_layer[i] == 0) score -= 1.;
    if(nhits_layer[i] > 0) score += 1./(nhits_layer[i]);
  }
  score = (score/N_HCAL_LAYERS + 1.)/2.;
  return score;
}

float MIP_Likeness_total(float nhits_e_layer[N_ECAL_LAYERS], float nhits_h_layer[N_HCAL_LAYERS]) {
  float score = 0.;
  float score_e = 0.;
  float score_h = 0.;
  for(int i=0; i<N_ECAL_LAYERS; i++) {
    if(nhits_e_layer[i] == 0) score_e -= 1.;
    if(nhits_e_layer[i] > 0) score_e += 1./(nhits_e_layer[i]);
  }
  for(int i=0; i<N_HCAL_LAYERS; i++) {
    if(nhits_h_layer[i] == 0) score_h -= 1.;
    if(nhits_h_layer[i] > 0) score_h += 1./(nhits_h_layer[i]);
  }
  score = (score_e/N_ECAL_LAYERS + score_h/N_HCAL_LAYERS + 2.)/4.;
  return score;
}

void hcal_ecal_ratio(float &nhit_ratio, float &sume_ratio, float &weight_ratio, int nhit_e, float sume_e, float weight_e, int nhit_h, float sume_h, float weight_h) {
  //I'll sum 1 hit in each detector to escape problems with 0, simple.
  nhit_ratio = (nhit_h+1.)/(nhit_e+1.);
  sume_ratio = (sume_h+1.)/(sume_e+1.);
  weight_ratio = (weight_h+1.)/(weight_e+1.);
  
  // Cut down
  //if(nhit_ratio>100.) nhit_ratio = 100.;
  //if(sume_ratio>100.) sume_ratio = 100.;
  //if(weight_ratio>100.) weight_ratio = 100.;

  return 0;
}

void is_interaction(float &ecal_int, float &hcal_int, float &total_int, int nhit_e, int nhit_h) {
  ecal_int = 0.;
  hcal_int = 0.;
  total_int = 0.;
  if(nhit_e > 0) ecal_int = 1.;
  if(nhit_h > 0) hcal_int = 1.;
  if((nhit_e > 0) or (nhit_h > 0)) total_int = 1.;
  return 0;
}

struct hitpair
{
  float hit_rs;
  float hit_es;
};

bool CompareHitsR(const hitpair hit1, const hitpair hit2) { return hit1.hit_rs < hit2.hit_rs; }

float moliere_ecal(vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses,
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

float moliere_hcal(vector<float> * hit_energy, vector<int> *hit_slab, TVectorD S_thicknesses,
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
      weighte += hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
      wx += hit_x->at(j) * hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
      wy += hit_y->at(j) * hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
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
      hits_vector.push_back({static_cast<float>(r),static_cast<float>(hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618)});
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

float hits_max_distance(vector<int> *hit_slab, vector<float> * hit_x, vector<float> * hit_y,
			     vector<int> * hit_isMasked, bool masked=false) {
  float max_distance = 0.;
  //cout<<"DEBUG1"<<endl;
  //cout<<"hit_slab.size(): "<<hit_slab->size()<<endl;
  if(hit_slab->size() > 1){
    for (int i = 0; i < hit_slab->size(); i++){
      if (masked && hit_isMasked->at(i) == 1) continue;
      //  cout<<"DEBUG2"<<endl;
      int layer1 = hit_slab->at(i);
      float hit1_x = hit_x->at(i);
      float hit1_y = hit_y->at(i);
      //cout<<"DEBUG2.1"<<endl;
      for (int j = i+1; j < hit_slab->size()-1; j++){
	//cout<<"DEBUG3"<<endl;
	//if(j==i) continue;
	if (masked && hit_isMasked->at(j) == 1) continue;
	int layer2 = hit_slab->at(j);
	//cout<<"DEBUG3.1"<<endl;
	if(layer2 != layer1) continue;
	float hit2_x = hit_x->at(j);
	float hit2_y = hit_y->at(j);
	float distance = pow(pow(hit2_x-hit1_x,2)+pow(hit2_y-hit1_y,2),0.5);
	if(distance>max_distance) max_distance = distance;
      }
    }
  }
  //  cout<<"DEBUG4"<<endl;
  return max_distance;
}


float moliere_total(vector<float> * e_hit_energy, vector<int> *e_hit_slab, TVectorD W_thicknesses,
                   vector<float> * e_hit_x, vector<float> * e_hit_y, vector<float> * e_hit_z, vector<int> * e_hit_isMasked,
		   vector<float> * h_hit_energy, vector<int> *h_hit_slab, TVectorD S_thicknesses,
                   vector<float> * h_hit_x, vector<float> * h_hit_y, vector<float> * h_hit_z, vector<int> * h_hit_isMasked, 
		   bool masked=false, float containment = 0.90, bool is_shower=false) {
  float mol_rad = 0.;
  float weighte = 0.;
  float wx = 0.; float wy = 0.; float wz = 0.;
  
  if(is_shower==true){
      vector<hitpair> hits_vector;
      if(e_hit_energy->size() > 0){
	for (int j = 0; j < e_hit_energy->size(); j++){
	  if (masked && e_hit_isMasked->at(j) == 1) continue;
	  weighte += e_hit_energy->at(j) * S_thicknesses[e_hit_slab->at(j)]/3.5;
	  wx += e_hit_x->at(j) * e_hit_energy->at(j) * S_thicknesses[e_hit_slab->at(j)]/3.5;
	  wy += e_hit_y->at(j) * e_hit_energy->at(j) * S_thicknesses[e_hit_slab->at(j)]/3.5;
	}
      }
      
      if(h_hit_energy->size() > 0){
	for (int j = 0; j < h_hit_energy->size(); j++){
	  if (masked && h_hit_isMasked->at(j) == 1) continue;
	  weighte += h_hit_energy->at(j) * S_thicknesses[h_hit_slab->at(j)]/17.618;
	  wx += h_hit_x->at(j) * h_hit_energy->at(j) * S_thicknesses[h_hit_slab->at(j)]/17.618;
	  wy += h_hit_y->at(j) * h_hit_energy->at(j) * S_thicknesses[h_hit_slab->at(j)]/17.618;
	}
      }
      
    float bary_x = wx / weighte;
    float bary_y = wy / weighte;
    if( (weighte == 0) or (weighte < 0) ){
      bary_x = 0.;
      bary_y = 0.;
    }

    if(e_hit_energy->size() > 0){
      for (int j = 0; j < e_hit_energy->size(); j++){
	if (masked && e_hit_isMasked->at(j) == 1) continue;
	float r = pow(pow((e_hit_x->at(j) - bary_x) , 2) + pow((e_hit_y->at(j) - bary_y), 2), 0.5);
	hits_vector.push_back({static_cast<float>(r),static_cast<float>(e_hit_energy->at(j) * W_thicknesses[e_hit_slab->at(j)]/3.5)});
      }
    }
    
    if(h_hit_energy->size() > 0){
      for (int j = 0; j < h_hit_energy->size(); j++){
        if (masked && h_hit_isMasked->at(j) == 1) continue;
        float r = pow(pow((h_hit_x->at(j) - bary_x) , 2) + pow((h_hit_y->at(j) - bary_y), 2), 0.5);
        hits_vector.push_back({static_cast<float>(r),static_cast<float>(h_hit_energy->at(j) * S_thicknesses[h_hit_slab->at(j)]/17.618)});
      }
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

void radius_layer_ecal(float mol_per_layer[N_ECAL_LAYERS], vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses,
                  vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z,
                  vector<int> * hit_isMasked, bool masked=false, float containment = 0.90, bool is_shower=false) {

  if((hit_energy->size() > 0) and (is_shower==true)){

    for(int ilayer = 0; ilayer < N_ECAL_LAYERS; ilayer++){
      float mol_rad = 0.;
      float weighte = 0.;
      float wx = 0.; float wy = 0.;
      float r = 0.;
      int nhit_layer = 0;
  
      vector<hitpair> hits_vector;
      
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
	std::sort(hits_vector.begin(), hits_vector.end(), CompareHitsR);

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

void radius_layer_hcal(float mol_per_layer[N_HCAL_LAYERS], vector<float> * hit_energy, vector<int> *hit_slab, TVectorD S_thicknesses,
		       vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z,
		       vector<int> * hit_isMasked, bool masked=false, float containment = 0.90, bool is_shower=false) {

  if((hit_energy->size() > 0) and (is_shower==true)){

    for(int ilayer = 0; ilayer < N_HCAL_LAYERS; ilayer++){
      float mol_rad = 0.;
      float weighte = 0.;
      float wx = 0.; float wy = 0.;
      float r = 0.;
      int nhit_layer = 0;

      vector<hitpair> hits_vector;

      for (int j = 0; j < hit_energy->size(); j++){
        if (masked && hit_isMasked->at(j) == 1) continue;
        if(hit_slab->at(j) == ilayer) {
          weighte += hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
          wx += hit_x->at(j) * hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
          wy += hit_y->at(j) * hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
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
          hits_vector.push_back({static_cast<float>(r),static_cast<float>(hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618)});
        }
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
        if(mol_i<0) mol_i=0.;
        mol_per_layer[ilayer] = hits_vector.at(mol_i).hit_rs;
      }
      if(nhit_layer < 3) mol_per_layer[ilayer] = 0.;
      if(mol_per_layer[ilayer] < 1) mol_per_layer[ilayer] = 0.;
      hits_vector.clear();
    }
  }
  return 0;
}

void barycenter_ecal(vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses,
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

void barycenter_hcal(vector<float> * hit_energy, vector<int> *hit_slab, TVectorD S_thicknesses,
		     vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z, float bar_xyzr[4],
		     vector<int> * hit_isMasked, bool masked=false, bool is_shower=false) {

  float sume = 0.;
  float wx = 0.; float wy = 0.; float wz = 0.;
  float bary_x = 0., bary_y = 0., bary_z = 0.;
  if(hit_energy->size() > 0){
    for (int j = 0; j < hit_energy->size(); j++){
      if (masked && hit_isMasked->at(j) == 1) continue;
      sume += hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
      wx += hit_x->at(j) * hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
      wy += hit_y->at(j) * hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
      wz += hit_z->at(j) * hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
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

void barycenter_total(vector<float> * e_hit_energy, vector<int> *e_hit_slab, TVectorD W_thicknesses,
		      vector<float> * e_hit_x, vector<float> * e_hit_y, vector<float> * e_hit_z, vector<int> * e_hit_isMasked,
		      vector<float> * h_hit_energy, vector<int> *h_hit_slab, TVectorD S_thicknesses, 
		      vector<float> * h_hit_x, vector<float> * h_hit_y, vector<float> * h_hit_z, vector<int> * h_hit_isMasked,
		      float bar_xyzr[4], bool masked=false) {

  float sume = 0.;
  float wx = 0.; float wy = 0.; float wz = 0.;
  float bary_x = 0., bary_y = 0., bary_z = 0.;
  if(e_hit_energy->size() > 0){
    for (int j = 0; j < e_hit_energy->size(); j++){
      if (masked && e_hit_isMasked->at(j) == 1) continue;
      sume += e_hit_energy->at(j) * W_thicknesses[e_hit_slab->at(j)]/3.5;
      wx += e_hit_x->at(j) * e_hit_energy->at(j) * W_thicknesses[e_hit_slab->at(j)]/3.5;
      wy += e_hit_y->at(j) * e_hit_energy->at(j) * W_thicknesses[e_hit_slab->at(j)]/3.5;
      wz += e_hit_z->at(j) * e_hit_energy->at(j) * W_thicknesses[e_hit_slab->at(j)]/3.5;
    }
  }
  if(h_hit_energy->size() > 0){
    for (int j = 0; j < h_hit_energy->size(); j++){
      if (masked && h_hit_isMasked->at(j) == 1) continue;
      sume += h_hit_energy->at(j) * S_thicknesses[h_hit_slab->at(j)]/17.618;
      wx += h_hit_x->at(j) * h_hit_energy->at(j) * S_thicknesses[h_hit_slab->at(j)]/17.618;
      wy += h_hit_y->at(j) * h_hit_energy->at(j) * S_thicknesses[h_hit_slab->at(j)]/17.618;
      wz += h_hit_z->at(j) * h_hit_energy->at(j) * S_thicknesses[h_hit_slab->at(j)]/17.618;
    }
  }
  if((e_hit_energy->size() > 0) or (h_hit_energy->size() > 0)){
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

void bary_layer_ecal(float blv[N_ECAL_LAYERS][3], vector<float> * hit_energy, vector<int> *hit_slab,
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

void bary_layer_hcal(float blv[N_HCAL_LAYERS][3], vector<float> * hit_energy, vector<int> *hit_slab,
		     TVectorD S_thicknesses, vector<float> * hit_x, vector<float> * hit_y,
		     vector<float> * hit_z, vector<int> * hit_isMasked, bool masked=false, bool is_shower=false) {

  float sume_w[N_HCAL_LAYERS];
  float wx[N_HCAL_LAYERS]; float wy[N_HCAL_LAYERS];
  float bary_x[N_HCAL_LAYERS]; float bary_y[N_HCAL_LAYERS];

  for (int ilayer=0; ilayer < N_HCAL_LAYERS; ilayer++){
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
    for (int ilayer=0; ilayer < N_HCAL_LAYERS; ilayer++) {
      if(hit_slab->size() > 0){
        for (int j = 0; j < hit_energy->size(); j++) {
          if (masked && hit_isMasked->at(j) == 1) continue;
          if (hit_slab->at(j) == 0) {
            sume_w[ilayer] += hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
            wx[ilayer] += hit_x->at(j) * hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
            wy[ilayer] += hit_y->at(j) * hit_energy->at(j) * S_thicknesses[hit_slab->at(j)]/17.618;
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
    
    double test_e[N_ENERGIES]={6.};
    TVectorD energies(N_ENERGIES, test_e);
    TVectorD energies_tr(N_ENERGIES);
    for (int j = 0; j < N_ENERGIES; j++) energies_tr[j] = 1/TMath::Sqrt(energies[j]);
    // For conf6
    double W[N_ECAL_LAYERS] = {4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6}; 
    TVectorD W_thicknesses(N_ECAL_LAYERS, W);
    // Steel thickness preshower: 0.5mm
    // Steel thickness per layer (next 38): (0.5 + 17.2 + 0.5) = 18.2mm
    // Steel thickness before bad layer: (0.5 + 17.2 + 17.2 + 0.5) = 35.4mm
    // Steel thickness before tokyo layer: (0.5 + 17.2 + 17.2 + 0.5) = 35.4mm

    double S[N_HCAL_LAYERS];
    S[0] = 0.5;
    for(int i = 1; i<39; i++) {
      S[i] = 0.5 + 17.2 + 0.5;
    }
    S[39] = 0.5 + 17.2 + 17.2 + 0.5;
    S[40] = 0.5 + 17.2 + 17.2 + 0.5;
    TVectorD S_thicknesses(N_HCAL_LAYERS, S);

    string filenames_ecal[N_ENERGIES];
    string filenames_hcal[N_ENERGIES];
    string base_path = "/lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/combined/submit_jobs/combined_LCIO2build_folder/combined_LCIO2build_output/";
    for (int j = 0; j < N_ENERGIES; j++){
      std::string str = std::to_string (energies[j]);
      str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
      str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
      cout<<str<<endl;
      filenames_ecal[j] = base_path + "output_LCIO2Build_TB2022-06_"+particle+"_" + str +  "GeV.root";
      filenames_hcal[j] = base_path + "output_AHCALLCIO2Build_TB2022-06_"+particle+"_" + str +  "GeV.root";
    }
    TString result_name = "resolution_"+particle+"_result.root" ;
    
    TFile f(result_name, "recreate");
    TTree *outtree = new TTree("ntp","NTuples");

    // branches definitions
    // We will save both the branches and the histos
    Float_t b_ecal_interaction;
    Float_t b_ecal_nhit, b_ecal_sume, b_ecal_weighte, b_ecal_bar_x, b_ecal_bar_y, b_ecal_bar_z, b_ecal_bar_r;
    Float_t b_ecal_mol;
    Float_t b_ecal_MIP_Likeness;
    Float_t b_ecal_hits_max_distance;
    Float_t b_ecal_radius90_layer_0, b_ecal_radius90_layer_1, b_ecal_radius90_layer_2, b_ecal_radius90_layer_3, b_ecal_radius90_layer_4, b_ecal_radius90_layer_5, b_ecal_radius90_layer_6, b_ecal_radius90_layer_7, b_ecal_radius90_layer_8, b_ecal_radius90_layer_9, b_ecal_radius90_layer_10, b_ecal_radius90_layer_11, b_ecal_radius90_layer_12, b_ecal_radius90_layer_13, b_ecal_radius90_layer_14;
    Float_t b_ecal_bar_x_layer_0, b_ecal_bar_x_layer_1, b_ecal_bar_x_layer_2, b_ecal_bar_x_layer_3, b_ecal_bar_x_layer_4, b_ecal_bar_x_layer_5, b_ecal_bar_x_layer_6, b_ecal_bar_x_layer_7, b_ecal_bar_x_layer_8, b_ecal_bar_x_layer_9, b_ecal_bar_x_layer_10, b_ecal_bar_x_layer_11, b_ecal_bar_x_layer_12, b_ecal_bar_x_layer_13, b_ecal_bar_x_layer_14;
    Float_t b_ecal_bar_y_layer_0, b_ecal_bar_y_layer_1, b_ecal_bar_y_layer_2, b_ecal_bar_y_layer_3, b_ecal_bar_y_layer_4, b_ecal_bar_y_layer_5, b_ecal_bar_y_layer_6, b_ecal_bar_y_layer_7, b_ecal_bar_y_layer_8, b_ecal_bar_y_layer_9, b_ecal_bar_y_layer_10, b_ecal_bar_y_layer_11, b_ecal_bar_y_layer_12, b_ecal_bar_y_layer_13, b_ecal_bar_y_layer_14;
    Float_t b_ecal_bar_r_layer_0, b_ecal_bar_r_layer_1, b_ecal_bar_r_layer_2, b_ecal_bar_r_layer_3, b_ecal_bar_r_layer_4, b_ecal_bar_r_layer_5, b_ecal_bar_r_layer_6, b_ecal_bar_r_layer_7, b_ecal_bar_r_layer_8, b_ecal_bar_r_layer_9, b_ecal_bar_r_layer_10, b_ecal_bar_r_layer_11, b_ecal_bar_r_layer_12, b_ecal_bar_r_layer_13, b_ecal_bar_r_layer_14;
    Float_t b_ecal_shower_nhit_max_layer, b_ecal_shower_nhit_start_layer, b_ecal_shower_nhit_end_layer, b_ecal_shower_nhit_start_10_layer, b_ecal_shower_nhit_end_10_layer, b_ecal_shower_nhit_average, b_ecal_shower_nhit_max;
    Float_t b_ecal_shower_sume_max_layer, b_ecal_shower_sume_start_layer, b_ecal_shower_sume_end_layer, b_ecal_shower_sume_start_10_layer, b_ecal_shower_sume_end_10_layer, b_ecal_shower_sume_average, b_ecal_shower_sume_max;
    Float_t b_ecal_shower_weighte_max_layer, b_ecal_shower_weighte_start_layer, b_ecal_shower_weighte_end_layer, b_ecal_shower_weighte_start_10_layer, b_ecal_shower_weighte_end_10_layer, b_ecal_shower_weighte_average, b_ecal_shower_weighte_max;
    Float_t b_ecal_nhit_layer_0, b_ecal_nhit_layer_1, b_ecal_nhit_layer_2, b_ecal_nhit_layer_3, b_ecal_nhit_layer_4, b_ecal_nhit_layer_5, b_ecal_nhit_layer_6, b_ecal_nhit_layer_7, b_ecal_nhit_layer_8, b_ecal_nhit_layer_9, b_ecal_nhit_layer_10, b_ecal_nhit_layer_11, b_ecal_nhit_layer_12, b_ecal_nhit_layer_13, b_ecal_nhit_layer_14;
    Float_t b_ecal_nhit_layer_n_0, b_ecal_nhit_layer_n_1, b_ecal_nhit_layer_n_2, b_ecal_nhit_layer_n_3, b_ecal_nhit_layer_n_4, b_ecal_nhit_layer_n_5, b_ecal_nhit_layer_n_6, b_ecal_nhit_layer_n_7, b_ecal_nhit_layer_n_8, b_ecal_nhit_layer_n_9, b_ecal_nhit_layer_n_10, b_ecal_nhit_layer_n_11, b_ecal_nhit_layer_n_12, b_ecal_nhit_layer_n_13, b_ecal_nhit_layer_n_14;
    Float_t b_ecal_weighte_layer_0, b_ecal_weighte_layer_1, b_ecal_weighte_layer_2, b_ecal_weighte_layer_3, b_ecal_weighte_layer_4, b_ecal_weighte_layer_5, b_ecal_weighte_layer_6, b_ecal_weighte_layer_7, b_ecal_weighte_layer_8, b_ecal_weighte_layer_9, b_ecal_weighte_layer_10, b_ecal_weighte_layer_11, b_ecal_weighte_layer_12, b_ecal_weighte_layer_13, b_ecal_weighte_layer_14;
    Float_t b_ecal_weighte_layer_n_0, b_ecal_weighte_layer_n_1, b_ecal_weighte_layer_n_2, b_ecal_weighte_layer_n_3, b_ecal_weighte_layer_n_4, b_ecal_weighte_layer_n_5, b_ecal_weighte_layer_n_6, b_ecal_weighte_layer_n_7, b_ecal_weighte_layer_n_8, b_ecal_weighte_layer_n_9, b_ecal_weighte_layer_n_10, b_ecal_weighte_layer_n_11, b_ecal_weighte_layer_n_12, b_ecal_weighte_layer_n_13, b_ecal_weighte_layer_n_14;    
    Float_t b_ecal_sume_layer_0, b_ecal_sume_layer_1, b_ecal_sume_layer_2, b_ecal_sume_layer_3, b_ecal_sume_layer_4, b_ecal_sume_layer_5, b_ecal_sume_layer_6, b_ecal_sume_layer_7, b_ecal_sume_layer_8, b_ecal_sume_layer_9, b_ecal_sume_layer_10, b_ecal_sume_layer_11, b_ecal_sume_layer_12, b_ecal_sume_layer_13, b_ecal_sume_layer_14;
    Float_t b_ecal_sume_layer_n_0, b_ecal_sume_layer_n_1, b_ecal_sume_layer_n_2, b_ecal_sume_layer_n_3, b_ecal_sume_layer_n_4, b_ecal_sume_layer_n_5, b_ecal_sume_layer_n_6, b_ecal_sume_layer_n_7, b_ecal_sume_layer_n_8, b_ecal_sume_layer_n_9, b_ecal_sume_layer_n_10, b_ecal_sume_layer_n_11, b_ecal_sume_layer_n_12, b_ecal_sume_layer_n_13, b_ecal_sume_layer_n_14;

    Float_t b_hcal_interaction;
    Float_t b_hcal_nhit, b_hcal_sume, b_hcal_weighte, b_hcal_bar_x, b_hcal_bar_y, b_hcal_bar_z, b_hcal_bar_r;
    Float_t b_hcal_mol;
    Float_t b_hcal_MIP_Likeness;
    Float_t b_hcal_hits_max_distance;
    Float_t b_hcal_radius90_layer_0, b_hcal_radius90_layer_1, b_hcal_radius90_layer_2, b_hcal_radius90_layer_3, b_hcal_radius90_layer_4, b_hcal_radius90_layer_5, b_hcal_radius90_layer_6, b_hcal_radius90_layer_7, b_hcal_radius90_layer_8, b_hcal_radius90_layer_9, b_hcal_radius90_layer_10, b_hcal_radius90_layer_11, b_hcal_radius90_layer_12, b_hcal_radius90_layer_13, b_hcal_radius90_layer_14;
    Float_t b_hcal_bar_x_layer_0, b_hcal_bar_x_layer_1, b_hcal_bar_x_layer_2, b_hcal_bar_x_layer_3, b_hcal_bar_x_layer_4, b_hcal_bar_x_layer_5, b_hcal_bar_x_layer_6, b_hcal_bar_x_layer_7, b_hcal_bar_x_layer_8, b_hcal_bar_x_layer_9, b_hcal_bar_x_layer_10, b_hcal_bar_x_layer_11, b_hcal_bar_x_layer_12, b_hcal_bar_x_layer_13, b_hcal_bar_x_layer_14;
    Float_t b_hcal_bar_y_layer_0, b_hcal_bar_y_layer_1, b_hcal_bar_y_layer_2, b_hcal_bar_y_layer_3, b_hcal_bar_y_layer_4, b_hcal_bar_y_layer_5, b_hcal_bar_y_layer_6, b_hcal_bar_y_layer_7, b_hcal_bar_y_layer_8, b_hcal_bar_y_layer_9, b_hcal_bar_y_layer_10, b_hcal_bar_y_layer_11, b_hcal_bar_y_layer_12, b_hcal_bar_y_layer_13, b_hcal_bar_y_layer_14;
    Float_t b_hcal_bar_r_layer_0, b_hcal_bar_r_layer_1, b_hcal_bar_r_layer_2, b_hcal_bar_r_layer_3, b_hcal_bar_r_layer_4, b_hcal_bar_r_layer_5, b_hcal_bar_r_layer_6, b_hcal_bar_r_layer_7, b_hcal_bar_r_layer_8, b_hcal_bar_r_layer_9, b_hcal_bar_r_layer_10, b_hcal_bar_r_layer_11, b_hcal_bar_r_layer_12, b_hcal_bar_r_layer_13, b_hcal_bar_r_layer_14;
    Float_t b_hcal_shower_nhit_max_layer, b_hcal_shower_nhit_start_layer, b_hcal_shower_nhit_end_layer, b_hcal_shower_nhit_start_10_layer, b_hcal_shower_nhit_end_10_layer, b_hcal_shower_nhit_average, b_hcal_shower_nhit_max;
    Float_t b_hcal_shower_sume_max_layer, b_hcal_shower_sume_start_layer, b_hcal_shower_sume_end_layer, b_hcal_shower_sume_start_10_layer, b_hcal_shower_sume_end_10_layer, b_hcal_shower_sume_average, b_hcal_shower_sume_max;
    Float_t b_hcal_shower_weighte_max_layer, b_hcal_shower_weighte_start_layer, b_hcal_shower_weighte_end_layer, b_hcal_shower_weighte_start_10_layer, b_hcal_shower_weighte_end_10_layer, b_hcal_shower_weighte_average, b_hcal_shower_weighte_max;

    Float_t b_hcal_nhit_layer_0, b_hcal_nhit_layer_1, b_hcal_nhit_layer_2, b_hcal_nhit_layer_3, b_hcal_nhit_layer_4, b_hcal_nhit_layer_5, b_hcal_nhit_layer_6, b_hcal_nhit_layer_7, b_hcal_nhit_layer_8, b_hcal_nhit_layer_9, b_hcal_nhit_layer_10, b_hcal_nhit_layer_11, b_hcal_nhit_layer_12, b_hcal_nhit_layer_13, b_hcal_nhit_layer_14;
    Float_t b_hcal_nhit_layer_15, b_hcal_nhit_layer_16, b_hcal_nhit_layer_17, b_hcal_nhit_layer_18, b_hcal_nhit_layer_19, b_hcal_nhit_layer_20, b_hcal_nhit_layer_21, b_hcal_nhit_layer_22, b_hcal_nhit_layer_23, b_hcal_nhit_layer_24, b_hcal_nhit_layer_25, b_hcal_nhit_layer_26, b_hcal_nhit_layer_27, b_hcal_nhit_layer_28, b_hcal_nhit_layer_29;
    Float_t b_hcal_nhit_layer_30, b_hcal_nhit_layer_31, b_hcal_nhit_layer_32, b_hcal_nhit_layer_33, b_hcal_nhit_layer_34, b_hcal_nhit_layer_35, b_hcal_nhit_layer_36, b_hcal_nhit_layer_37, b_hcal_nhit_layer_38, b_hcal_nhit_layer_39, b_hcal_nhit_layer_40;

    Float_t b_hcal_nhit_layer_n_0, b_hcal_nhit_layer_n_1, b_hcal_nhit_layer_n_2, b_hcal_nhit_layer_n_3, b_hcal_nhit_layer_n_4, b_hcal_nhit_layer_n_5, b_hcal_nhit_layer_n_6, b_hcal_nhit_layer_n_7, b_hcal_nhit_layer_n_8, b_hcal_nhit_layer_n_9, b_hcal_nhit_layer_n_10, b_hcal_nhit_layer_n_11, b_hcal_nhit_layer_n_12, b_hcal_nhit_layer_n_13, b_hcal_nhit_layer_n_14;
    Float_t b_hcal_nhit_layer_n_15, b_hcal_nhit_layer_n_16, b_hcal_nhit_layer_n_17, b_hcal_nhit_layer_n_18, b_hcal_nhit_layer_n_19, b_hcal_nhit_layer_n_20, b_hcal_nhit_layer_n_21, b_hcal_nhit_layer_n_22, b_hcal_nhit_layer_n_23, b_hcal_nhit_layer_n_24, b_hcal_nhit_layer_n_25, b_hcal_nhit_layer_n_26, b_hcal_nhit_layer_n_27, b_hcal_nhit_layer_n_28, b_hcal_nhit_layer_n_29;
    Float_t b_hcal_nhit_layer_n_30, b_hcal_nhit_layer_n_31, b_hcal_nhit_layer_n_32, b_hcal_nhit_layer_n_33, b_hcal_nhit_layer_n_34, b_hcal_nhit_layer_n_35, b_hcal_nhit_layer_n_36, b_hcal_nhit_layer_n_37, b_hcal_nhit_layer_n_38, b_hcal_nhit_layer_n_39, b_hcal_nhit_layer_n_40;

    Float_t b_hcal_weighte_layer_0, b_hcal_weighte_layer_1, b_hcal_weighte_layer_2, b_hcal_weighte_layer_3, b_hcal_weighte_layer_4, b_hcal_weighte_layer_5, b_hcal_weighte_layer_6, b_hcal_weighte_layer_7, b_hcal_weighte_layer_8, b_hcal_weighte_layer_9, b_hcal_weighte_layer_10, b_hcal_weighte_layer_11, b_hcal_weighte_layer_12, b_hcal_weighte_layer_13, b_hcal_weighte_layer_14;
    Float_t b_hcal_weighte_layer_15, b_hcal_weighte_layer_16, b_hcal_weighte_layer_17, b_hcal_weighte_layer_18, b_hcal_weighte_layer_19, b_hcal_weighte_layer_20, b_hcal_weighte_layer_21, b_hcal_weighte_layer_22, b_hcal_weighte_layer_23, b_hcal_weighte_layer_24, b_hcal_weighte_layer_25, b_hcal_weighte_layer_26, b_hcal_weighte_layer_27, b_hcal_weighte_layer_28, b_hcal_weighte_layer_29;
    Float_t b_hcal_weighte_layer_30, b_hcal_weighte_layer_31, b_hcal_weighte_layer_32, b_hcal_weighte_layer_33, b_hcal_weighte_layer_34, b_hcal_weighte_layer_35, b_hcal_weighte_layer_36, b_hcal_weighte_layer_37, b_hcal_weighte_layer_38, b_hcal_weighte_layer_39, b_hcal_weighte_layer_40;

    Float_t b_hcal_weighte_layer_n_0, b_hcal_weighte_layer_n_1, b_hcal_weighte_layer_n_2, b_hcal_weighte_layer_n_3, b_hcal_weighte_layer_n_4, b_hcal_weighte_layer_n_5, b_hcal_weighte_layer_n_6, b_hcal_weighte_layer_n_7, b_hcal_weighte_layer_n_8, b_hcal_weighte_layer_n_9, b_hcal_weighte_layer_n_10, b_hcal_weighte_layer_n_11, b_hcal_weighte_layer_n_12, b_hcal_weighte_layer_n_13, b_hcal_weighte_layer_n_14;
    Float_t b_hcal_weighte_layer_n_15, b_hcal_weighte_layer_n_16, b_hcal_weighte_layer_n_17, b_hcal_weighte_layer_n_18, b_hcal_weighte_layer_n_19, b_hcal_weighte_layer_n_20, b_hcal_weighte_layer_n_21, b_hcal_weighte_layer_n_22, b_hcal_weighte_layer_n_23, b_hcal_weighte_layer_n_24, b_hcal_weighte_layer_n_25, b_hcal_weighte_layer_n_26, b_hcal_weighte_layer_n_27, b_hcal_weighte_layer_n_28, b_hcal_weighte_layer_n_29;
    Float_t b_hcal_weighte_layer_n_30, b_hcal_weighte_layer_n_31, b_hcal_weighte_layer_n_32, b_hcal_weighte_layer_n_33, b_hcal_weighte_layer_n_34, b_hcal_weighte_layer_n_35, b_hcal_weighte_layer_n_36, b_hcal_weighte_layer_n_37, b_hcal_weighte_layer_n_38, b_hcal_weighte_layer_n_39, b_hcal_weighte_layer_n_40;

    Float_t b_hcal_sume_layer_0, b_hcal_sume_layer_1, b_hcal_sume_layer_2, b_hcal_sume_layer_3, b_hcal_sume_layer_4, b_hcal_sume_layer_5, b_hcal_sume_layer_6, b_hcal_sume_layer_7, b_hcal_sume_layer_8, b_hcal_sume_layer_9, b_hcal_sume_layer_10, b_hcal_sume_layer_11, b_hcal_sume_layer_12, b_hcal_sume_layer_13, b_hcal_sume_layer_14;
    Float_t b_hcal_sume_layer_15, b_hcal_sume_layer_16, b_hcal_sume_layer_17, b_hcal_sume_layer_18, b_hcal_sume_layer_19, b_hcal_sume_layer_20, b_hcal_sume_layer_21, b_hcal_sume_layer_22, b_hcal_sume_layer_23, b_hcal_sume_layer_24, b_hcal_sume_layer_25, b_hcal_sume_layer_26, b_hcal_sume_layer_27, b_hcal_sume_layer_28, b_hcal_sume_layer_29;
    Float_t b_hcal_sume_layer_30, b_hcal_sume_layer_31, b_hcal_sume_layer_32, b_hcal_sume_layer_33, b_hcal_sume_layer_34, b_hcal_sume_layer_35, b_hcal_sume_layer_36, b_hcal_sume_layer_37, b_hcal_sume_layer_38, b_hcal_sume_layer_39, b_hcal_sume_layer_40;

    Float_t b_hcal_sume_layer_n_0, b_hcal_sume_layer_n_1, b_hcal_sume_layer_n_2, b_hcal_sume_layer_n_3, b_hcal_sume_layer_n_4, b_hcal_sume_layer_n_5, b_hcal_sume_layer_n_6, b_hcal_sume_layer_n_7, b_hcal_sume_layer_n_8, b_hcal_sume_layer_n_9, b_hcal_sume_layer_n_10, b_hcal_sume_layer_n_11, b_hcal_sume_layer_n_12, b_hcal_sume_layer_n_13, b_hcal_sume_layer_n_14;
    Float_t b_hcal_sume_layer_n_15, b_hcal_sume_layer_n_16, b_hcal_sume_layer_n_17, b_hcal_sume_layer_n_18, b_hcal_sume_layer_n_19, b_hcal_sume_layer_n_20, b_hcal_sume_layer_n_21, b_hcal_sume_layer_n_22, b_hcal_sume_layer_n_23, b_hcal_sume_layer_n_24, b_hcal_sume_layer_n_25, b_hcal_sume_layer_n_26, b_hcal_sume_layer_n_27, b_hcal_sume_layer_n_28, b_hcal_sume_layer_n_29;
    Float_t b_hcal_sume_layer_n_30, b_hcal_sume_layer_n_31, b_hcal_sume_layer_n_32, b_hcal_sume_layer_n_33, b_hcal_sume_layer_n_34, b_hcal_sume_layer_n_35, b_hcal_sume_layer_n_36, b_hcal_sume_layer_n_37, b_hcal_sume_layer_n_38, b_hcal_sume_layer_n_39, b_hcal_sume_layer_n_40;

    Float_t b_total_interaction;
    Float_t b_total_nhit, b_total_sume, b_total_weighte, b_total_bar_x, b_total_bar_y, b_total_bar_z, b_total_bar_r;
    Float_t b_total_mol;
    Float_t b_total_MIP_Likeness;
    Float_t b_total_hits_max_distance;
    Float_t b_total_radius90_layer_0, b_total_radius90_layer_1, b_total_radius90_layer_2, b_total_radius90_layer_3, b_total_radius90_layer_4, b_total_radius90_layer_5, b_total_radius90_layer_6, b_total_radius90_layer_7, b_total_radius90_layer_8, b_total_radius90_layer_9, b_total_radius90_layer_10, b_total_radius90_layer_11, b_total_radius90_layer_12, b_total_radius90_layer_13, b_total_radius90_layer_14;
    Float_t b_total_bar_x_layer_0, b_total_bar_x_layer_1, b_total_bar_x_layer_2, b_total_bar_x_layer_3, b_total_bar_x_layer_4, b_total_bar_x_layer_5, b_total_bar_x_layer_6, b_total_bar_x_layer_7, b_total_bar_x_layer_8, b_total_bar_x_layer_9, b_total_bar_x_layer_10, b_total_bar_x_layer_11, b_total_bar_x_layer_12, b_total_bar_x_layer_13, b_total_bar_x_layer_14;
    Float_t b_total_bar_y_layer_0, b_total_bar_y_layer_1, b_total_bar_y_layer_2, b_total_bar_y_layer_3, b_total_bar_y_layer_4, b_total_bar_y_layer_5, b_total_bar_y_layer_6, b_total_bar_y_layer_7, b_total_bar_y_layer_8, b_total_bar_y_layer_9, b_total_bar_y_layer_10, b_total_bar_y_layer_11, b_total_bar_y_layer_12, b_total_bar_y_layer_13, b_total_bar_y_layer_14;
    Float_t b_total_bar_r_layer_0, b_total_bar_r_layer_1, b_total_bar_r_layer_2, b_total_bar_r_layer_3, b_total_bar_r_layer_4, b_total_bar_r_layer_5, b_total_bar_r_layer_6, b_total_bar_r_layer_7, b_total_bar_r_layer_8, b_total_bar_r_layer_9, b_total_bar_r_layer_10, b_total_bar_r_layer_11, b_total_bar_r_layer_12, b_total_bar_r_layer_13, b_total_bar_r_layer_14;
    Float_t b_total_shower_nhit_max_layer, b_total_shower_nhit_start_layer, b_total_shower_nhit_end_layer, b_total_shower_nhit_start_10_layer, b_total_shower_nhit_end_10_layer, b_total_shower_nhit_average, b_total_shower_nhit_max;
    Float_t b_total_shower_sume_max_layer, b_total_shower_sume_start_layer, b_total_shower_sume_end_layer, b_total_shower_sume_start_10_layer, b_total_shower_sume_end_10_layer, b_total_shower_sume_average, b_total_shower_sume_max;
    Float_t b_total_shower_weighte_max_layer, b_total_shower_weighte_start_layer, b_total_shower_weighte_end_layer, b_total_shower_weighte_start_10_layer, b_total_shower_weighte_end_10_layer, b_total_shower_weighte_average, b_total_shower_weighte_max;
    Float_t b_total_nhit_layer_0, b_total_nhit_layer_1, b_total_nhit_layer_2, b_total_nhit_layer_3, b_total_nhit_layer_4, b_total_nhit_layer_5, b_total_nhit_layer_6, b_total_nhit_layer_7, b_total_nhit_layer_8, b_total_nhit_layer_9, b_total_nhit_layer_10, b_total_nhit_layer_11, b_total_nhit_layer_12, b_total_nhit_layer_13, b_total_nhit_layer_14;
    Float_t b_total_nhit_layer_15, b_total_nhit_layer_16, b_total_nhit_layer_17, b_total_nhit_layer_18, b_total_nhit_layer_19, b_total_nhit_layer_20, b_total_nhit_layer_21, b_total_nhit_layer_22, b_total_nhit_layer_23, b_total_nhit_layer_24, b_total_nhit_layer_25, b_total_nhit_layer_26, b_total_nhit_layer_27, b_total_nhit_layer_28, b_total_nhit_layer_29;
    Float_t b_total_nhit_layer_30, b_total_nhit_layer_31, b_total_nhit_layer_32, b_total_nhit_layer_33, b_total_nhit_layer_34, b_total_nhit_layer_35, b_total_nhit_layer_36, b_total_nhit_layer_37, b_total_nhit_layer_38, b_total_nhit_layer_39, b_total_nhit_layer_40, b_total_nhit_layer_41, b_total_nhit_layer_42, b_total_nhit_layer_43, b_total_nhit_layer_44;
    Float_t b_total_nhit_layer_45, b_total_nhit_layer_46, b_total_nhit_layer_47, b_total_nhit_layer_48, b_total_nhit_layer_49, b_total_nhit_layer_50, b_total_nhit_layer_51, b_total_nhit_layer_52, b_total_nhit_layer_53, b_total_nhit_layer_54, b_total_nhit_layer_55;

    Float_t b_total_nhit_layer_n_0, b_total_nhit_layer_n_1, b_total_nhit_layer_n_2, b_total_nhit_layer_n_3, b_total_nhit_layer_n_4, b_total_nhit_layer_n_5, b_total_nhit_layer_n_6, b_total_nhit_layer_n_7, b_total_nhit_layer_n_8, b_total_nhit_layer_n_9, b_total_nhit_layer_n_10, b_total_nhit_layer_n_11, b_total_nhit_layer_n_12, b_total_nhit_layer_n_13, b_total_nhit_layer_n_14;
    Float_t b_total_nhit_layer_n_15, b_total_nhit_layer_n_16, b_total_nhit_layer_n_17, b_total_nhit_layer_n_18, b_total_nhit_layer_n_19, b_total_nhit_layer_n_20, b_total_nhit_layer_n_21, b_total_nhit_layer_n_22, b_total_nhit_layer_n_23, b_total_nhit_layer_n_24, b_total_nhit_layer_n_25, b_total_nhit_layer_n_26, b_total_nhit_layer_n_27, b_total_nhit_layer_n_28, b_total_nhit_layer_n_29;
    Float_t b_total_nhit_layer_n_30, b_total_nhit_layer_n_31, b_total_nhit_layer_n_32, b_total_nhit_layer_n_33, b_total_nhit_layer_n_34, b_total_nhit_layer_n_35, b_total_nhit_layer_n_36, b_total_nhit_layer_n_37, b_total_nhit_layer_n_38, b_total_nhit_layer_n_39, b_total_nhit_layer_n_40, b_total_nhit_layer_n_41, b_total_nhit_layer_n_42, b_total_nhit_layer_n_43, b_total_nhit_layer_n_44;
    Float_t b_total_nhit_layer_n_45, b_total_nhit_layer_n_46, b_total_nhit_layer_n_47, b_total_nhit_layer_n_48, b_total_nhit_layer_n_49, b_total_nhit_layer_n_50, b_total_nhit_layer_n_51, b_total_nhit_layer_n_52, b_total_nhit_layer_n_53, b_total_nhit_layer_n_54, b_total_nhit_layer_n_55;

    Float_t b_total_weighte_layer_0, b_total_weighte_layer_1, b_total_weighte_layer_2, b_total_weighte_layer_3, b_total_weighte_layer_4, b_total_weighte_layer_5, b_total_weighte_layer_6, b_total_weighte_layer_7, b_total_weighte_layer_8, b_total_weighte_layer_9, b_total_weighte_layer_10, b_total_weighte_layer_11, b_total_weighte_layer_12, b_total_weighte_layer_13, b_total_weighte_layer_14;
    Float_t b_total_weighte_layer_15, b_total_weighte_layer_16, b_total_weighte_layer_17, b_total_weighte_layer_18, b_total_weighte_layer_19, b_total_weighte_layer_20, b_total_weighte_layer_21, b_total_weighte_layer_22, b_total_weighte_layer_23, b_total_weighte_layer_24, b_total_weighte_layer_25, b_total_weighte_layer_26, b_total_weighte_layer_27, b_total_weighte_layer_28, b_total_weighte_layer_29;
    Float_t b_total_weighte_layer_30, b_total_weighte_layer_31, b_total_weighte_layer_32, b_total_weighte_layer_33, b_total_weighte_layer_34, b_total_weighte_layer_35, b_total_weighte_layer_36, b_total_weighte_layer_37, b_total_weighte_layer_38, b_total_weighte_layer_39, b_total_weighte_layer_40, b_total_weighte_layer_41, b_total_weighte_layer_42, b_total_weighte_layer_43, b_total_weighte_layer_44;
    Float_t b_total_weighte_layer_45, b_total_weighte_layer_46, b_total_weighte_layer_47, b_total_weighte_layer_48, b_total_weighte_layer_49, b_total_weighte_layer_50, b_total_weighte_layer_51, b_total_weighte_layer_52, b_total_weighte_layer_53, b_total_weighte_layer_54, b_total_weighte_layer_55;

    Float_t b_total_weighte_layer_n_0, b_total_weighte_layer_n_1, b_total_weighte_layer_n_2, b_total_weighte_layer_n_3, b_total_weighte_layer_n_4, b_total_weighte_layer_n_5, b_total_weighte_layer_n_6, b_total_weighte_layer_n_7, b_total_weighte_layer_n_8, b_total_weighte_layer_n_9, b_total_weighte_layer_n_10, b_total_weighte_layer_n_11, b_total_weighte_layer_n_12, b_total_weighte_layer_n_13, b_total_weighte_layer_n_14;
    Float_t b_total_weighte_layer_n_15, b_total_weighte_layer_n_16, b_total_weighte_layer_n_17, b_total_weighte_layer_n_18, b_total_weighte_layer_n_19, b_total_weighte_layer_n_20, b_total_weighte_layer_n_21, b_total_weighte_layer_n_22, b_total_weighte_layer_n_23, b_total_weighte_layer_n_24, b_total_weighte_layer_n_25, b_total_weighte_layer_n_26, b_total_weighte_layer_n_27, b_total_weighte_layer_n_28, b_total_weighte_layer_n_29;
    Float_t b_total_weighte_layer_n_30, b_total_weighte_layer_n_31, b_total_weighte_layer_n_32, b_total_weighte_layer_n_33, b_total_weighte_layer_n_34, b_total_weighte_layer_n_35, b_total_weighte_layer_n_36, b_total_weighte_layer_n_37, b_total_weighte_layer_n_38, b_total_weighte_layer_n_39, b_total_weighte_layer_n_40, b_total_weighte_layer_n_41, b_total_weighte_layer_n_42, b_total_weighte_layer_n_43, b_total_weighte_layer_n_44;
    Float_t b_total_weighte_layer_n_45, b_total_weighte_layer_n_46, b_total_weighte_layer_n_47, b_total_weighte_layer_n_48, b_total_weighte_layer_n_49, b_total_weighte_layer_n_50, b_total_weighte_layer_n_51, b_total_weighte_layer_n_52, b_total_weighte_layer_n_53, b_total_weighte_layer_n_54, b_total_weighte_layer_n_55;

    Float_t b_total_sume_layer_0, b_total_sume_layer_1, b_total_sume_layer_2, b_total_sume_layer_3, b_total_sume_layer_4, b_total_sume_layer_5, b_total_sume_layer_6, b_total_sume_layer_7, b_total_sume_layer_8, b_total_sume_layer_9, b_total_sume_layer_10, b_total_sume_layer_11, b_total_sume_layer_12, b_total_sume_layer_13, b_total_sume_layer_14;
    Float_t b_total_sume_layer_15, b_total_sume_layer_16, b_total_sume_layer_17, b_total_sume_layer_18, b_total_sume_layer_19, b_total_sume_layer_20, b_total_sume_layer_21, b_total_sume_layer_22, b_total_sume_layer_23, b_total_sume_layer_24, b_total_sume_layer_25, b_total_sume_layer_26, b_total_sume_layer_27, b_total_sume_layer_28, b_total_sume_layer_29;
    Float_t b_total_sume_layer_30, b_total_sume_layer_31, b_total_sume_layer_32, b_total_sume_layer_33, b_total_sume_layer_34, b_total_sume_layer_35, b_total_sume_layer_36, b_total_sume_layer_37, b_total_sume_layer_38, b_total_sume_layer_39, b_total_sume_layer_40, b_total_sume_layer_41, b_total_sume_layer_42, b_total_sume_layer_43, b_total_sume_layer_44;
    Float_t b_total_sume_layer_45, b_total_sume_layer_46, b_total_sume_layer_47, b_total_sume_layer_48, b_total_sume_layer_49, b_total_sume_layer_50, b_total_sume_layer_51, b_total_sume_layer_52, b_total_sume_layer_53, b_total_sume_layer_54, b_total_sume_layer_55;

    Float_t b_total_sume_layer_n_0, b_total_sume_layer_n_1, b_total_sume_layer_n_2, b_total_sume_layer_n_3, b_total_sume_layer_n_4, b_total_sume_layer_n_5, b_total_sume_layer_n_6, b_total_sume_layer_n_7, b_total_sume_layer_n_8, b_total_sume_layer_n_9, b_total_sume_layer_n_10, b_total_sume_layer_n_11, b_total_sume_layer_n_12, b_total_sume_layer_n_13, b_total_sume_layer_n_14;
    Float_t b_total_sume_layer_n_15, b_total_sume_layer_n_16, b_total_sume_layer_n_17, b_total_sume_layer_n_18, b_total_sume_layer_n_19, b_total_sume_layer_n_20, b_total_sume_layer_n_21, b_total_sume_layer_n_22, b_total_sume_layer_n_23, b_total_sume_layer_n_24, b_total_sume_layer_n_25, b_total_sume_layer_n_26, b_total_sume_layer_n_27, b_total_sume_layer_n_28, b_total_sume_layer_n_29;
    Float_t b_total_sume_layer_n_30, b_total_sume_layer_n_31, b_total_sume_layer_n_32, b_total_sume_layer_n_33, b_total_sume_layer_n_34, b_total_sume_layer_n_35, b_total_sume_layer_n_36, b_total_sume_layer_n_37, b_total_sume_layer_n_38, b_total_sume_layer_n_39, b_total_sume_layer_n_40, b_total_sume_layer_n_41, b_total_sume_layer_n_42, b_total_sume_layer_n_43, b_total_sume_layer_n_44;
    Float_t b_total_sume_layer_n_45, b_total_sume_layer_n_46, b_total_sume_layer_n_47, b_total_sume_layer_n_48, b_total_sume_layer_n_49, b_total_sume_layer_n_50, b_total_sume_layer_n_51, b_total_sume_layer_n_52, b_total_sume_layer_n_53, b_total_sume_layer_n_54, b_total_sume_layer_n_55;

    Float_t b_nhit_ratio, b_sume_ratio, b_weighte_ratio;
        
    // Adding the branches
    outtree->Branch("ecal_interaction",&b_ecal_interaction,"b_ecal_interaction/F");
    outtree->Branch("ecal_nhit",&b_ecal_nhit,"b_ecal_nhit/F");
    outtree->Branch("ecal_sume",&b_ecal_sume,"b_ecal_sume/F");
    outtree->Branch("ecal_weighte",&b_ecal_weighte,"b_ecal_weighte/F");
    outtree->Branch("ecal_mol",&b_ecal_mol,"b_ecal_mol/F");
    outtree->Branch("ecal_MIP_Likeness",&b_ecal_MIP_Likeness,"b_ecal_MIP_Likeness/F");
    outtree->Branch("ecal_hits_max_distance",&b_ecal_hits_max_distance,"b_ecal_hits_max_distance/F");
    outtree->Branch("ecal_radius90_layer_0",&b_ecal_radius90_layer_0,"b_ecal_radius90_layer_0/F");
    outtree->Branch("ecal_radius90_layer_1",&b_ecal_radius90_layer_1,"b_ecal_radius90_layer_1/F");
    outtree->Branch("ecal_radius90_layer_2",&b_ecal_radius90_layer_2,"b_ecal_radius90_layer_2/F");
    outtree->Branch("ecal_radius90_layer_3",&b_ecal_radius90_layer_3,"b_ecal_radius90_layer_3/F");
    outtree->Branch("ecal_radius90_layer_4",&b_ecal_radius90_layer_4,"b_ecal_radius90_layer_4/F");
    outtree->Branch("ecal_radius90_layer_5",&b_ecal_radius90_layer_5,"b_ecal_radius90_layer_5/F");
    outtree->Branch("ecal_radius90_layer_6",&b_ecal_radius90_layer_6,"b_ecal_radius90_layer_6/F");
    outtree->Branch("ecal_radius90_layer_7",&b_ecal_radius90_layer_7,"b_ecal_radius90_layer_7/F");
    outtree->Branch("ecal_radius90_layer_8",&b_ecal_radius90_layer_8,"b_ecal_radius90_layer_8/F");
    outtree->Branch("ecal_radius90_layer_9",&b_ecal_radius90_layer_9,"b_ecal_radius90_layer_9/F");
    outtree->Branch("ecal_radius90_layer_10",&b_ecal_radius90_layer_10,"b_ecal_radius90_layer_10/F");
    outtree->Branch("ecal_radius90_layer_11",&b_ecal_radius90_layer_11,"b_ecal_radius90_layer_11/F");
    outtree->Branch("ecal_radius90_layer_12",&b_ecal_radius90_layer_12,"b_ecal_radius90_layer_12/F");
    outtree->Branch("ecal_radius90_layer_13",&b_ecal_radius90_layer_13,"b_ecal_radius90_layer_13/F");
    outtree->Branch("ecal_radius90_layer_14",&b_ecal_radius90_layer_14,"b_ecal_radius90_layer_14/F");
    outtree->Branch("ecal_bar_x",&b_ecal_bar_x,"b_ecal_bar_x/F");
    outtree->Branch("ecal_bar_y",&b_ecal_bar_y,"b_ecal_bar_y/F");
    outtree->Branch("ecal_bar_z",&b_ecal_bar_z,"b_ecal_bar_z/F");
    outtree->Branch("ecal_bar_r",&b_ecal_bar_r,"b_ecal_bar_r/F");
    outtree->Branch("ecal_bar_x_layer_0",&b_ecal_bar_x_layer_0,"b_ecal_bar_x_layer_0/F");
    outtree->Branch("ecal_bar_x_layer_1",&b_ecal_bar_x_layer_1,"b_ecal_bar_x_layer_1/F");
    outtree->Branch("ecal_bar_x_layer_2",&b_ecal_bar_x_layer_2,"b_ecal_bar_x_layer_2/F");
    outtree->Branch("ecal_bar_x_layer_3",&b_ecal_bar_x_layer_3,"b_ecal_bar_x_layer_3/F");
    outtree->Branch("ecal_bar_x_layer_4",&b_ecal_bar_x_layer_4,"b_ecal_bar_x_layer_4/F");
    outtree->Branch("ecal_bar_x_layer_5",&b_ecal_bar_x_layer_5,"b_ecal_bar_x_layer_5/F");
    outtree->Branch("ecal_bar_x_layer_6",&b_ecal_bar_x_layer_6,"b_ecal_bar_x_layer_6/F");
    outtree->Branch("ecal_bar_x_layer_7",&b_ecal_bar_x_layer_7,"b_ecal_bar_x_layer_7/F");
    outtree->Branch("ecal_bar_x_layer_8",&b_ecal_bar_x_layer_8,"b_ecal_bar_x_layer_8/F");
    outtree->Branch("ecal_bar_x_layer_9",&b_ecal_bar_x_layer_9,"b_ecal_bar_x_layer_9/F");
    outtree->Branch("ecal_bar_x_layer_10",&b_ecal_bar_x_layer_10,"b_ecal_bar_x_layer_10/F");
    outtree->Branch("ecal_bar_x_layer_11",&b_ecal_bar_x_layer_11,"b_ecal_bar_x_layer_11/F");
    outtree->Branch("ecal_bar_x_layer_12",&b_ecal_bar_x_layer_12,"b_ecal_bar_x_layer_12/F");
    outtree->Branch("ecal_bar_x_layer_13",&b_ecal_bar_x_layer_13,"b_ecal_bar_x_layer_13/F");
    outtree->Branch("ecal_bar_x_layer_14",&b_ecal_bar_x_layer_14,"b_ecal_bar_x_layer_14/F");
    outtree->Branch("ecal_bar_y_layer_0",&b_ecal_bar_y_layer_0,"b_ecal_bar_y_layer_0/F");
    outtree->Branch("ecal_bar_y_layer_1",&b_ecal_bar_y_layer_1,"b_ecal_bar_y_layer_1/F");
    outtree->Branch("ecal_bar_y_layer_2",&b_ecal_bar_y_layer_2,"b_ecal_bar_y_layer_2/F");
    outtree->Branch("ecal_bar_y_layer_3",&b_ecal_bar_y_layer_3,"b_ecal_bar_y_layer_3/F");
    outtree->Branch("ecal_bar_y_layer_4",&b_ecal_bar_y_layer_4,"b_ecal_bar_y_layer_4/F");
    outtree->Branch("ecal_bar_y_layer_5",&b_ecal_bar_y_layer_5,"b_ecal_bar_y_layer_5/F");
    outtree->Branch("ecal_bar_y_layer_6",&b_ecal_bar_y_layer_6,"b_ecal_bar_y_layer_6/F");
    outtree->Branch("ecal_bar_y_layer_7",&b_ecal_bar_y_layer_7,"b_ecal_bar_y_layer_7/F");
    outtree->Branch("ecal_bar_y_layer_8",&b_ecal_bar_y_layer_8,"b_ecal_bar_y_layer_8/F");
    outtree->Branch("ecal_bar_y_layer_9",&b_ecal_bar_y_layer_9,"b_ecal_bar_y_layer_9/F");
    outtree->Branch("ecal_bar_y_layer_10",&b_ecal_bar_y_layer_10,"b_ecal_bar_y_layer_10/F");
    outtree->Branch("ecal_bar_y_layer_11",&b_ecal_bar_y_layer_11,"b_ecal_bar_y_layer_11/F");
    outtree->Branch("ecal_bar_y_layer_12",&b_ecal_bar_y_layer_12,"b_ecal_bar_y_layer_12/F");
    outtree->Branch("ecal_bar_y_layer_13",&b_ecal_bar_y_layer_13,"b_ecal_bar_y_layer_13/F");
    outtree->Branch("ecal_bar_y_layer_14",&b_ecal_bar_y_layer_14,"b_ecal_bar_y_layer_14/F");
    outtree->Branch("ecal_bar_r_layer_0",&b_ecal_bar_r_layer_0,"b_ecal_bar_r_layer_0/F");
    outtree->Branch("ecal_bar_r_layer_1",&b_ecal_bar_r_layer_1,"b_ecal_bar_r_layer_1/F");
    outtree->Branch("ecal_bar_r_layer_2",&b_ecal_bar_r_layer_2,"b_ecal_bar_r_layer_2/F");
    outtree->Branch("ecal_bar_r_layer_3",&b_ecal_bar_r_layer_3,"b_ecal_bar_r_layer_3/F");
    outtree->Branch("ecal_bar_r_layer_4",&b_ecal_bar_r_layer_4,"b_ecal_bar_r_layer_4/F");
    outtree->Branch("ecal_bar_r_layer_5",&b_ecal_bar_r_layer_5,"b_ecal_bar_r_layer_5/F");
    outtree->Branch("ecal_bar_r_layer_6",&b_ecal_bar_r_layer_6,"b_ecal_bar_r_layer_6/F");
    outtree->Branch("ecal_bar_r_layer_7",&b_ecal_bar_r_layer_7,"b_ecal_bar_r_layer_7/F");
    outtree->Branch("ecal_bar_r_layer_8",&b_ecal_bar_r_layer_8,"b_ecal_bar_r_layer_8/F");
    outtree->Branch("ecal_bar_r_layer_9",&b_ecal_bar_r_layer_9,"b_ecal_bar_r_layer_9/F");
    outtree->Branch("ecal_bar_r_layer_10",&b_ecal_bar_r_layer_10,"b_ecal_bar_r_layer_10/F");
    outtree->Branch("ecal_bar_r_layer_11",&b_ecal_bar_r_layer_11,"b_ecal_bar_r_layer_11/F");
    outtree->Branch("ecal_bar_r_layer_12",&b_ecal_bar_r_layer_12,"b_ecal_bar_r_layer_12/F");
    outtree->Branch("ecal_bar_r_layer_13",&b_ecal_bar_r_layer_13,"b_ecal_bar_r_layer_13/F");
    outtree->Branch("ecal_bar_r_layer_14",&b_ecal_bar_r_layer_14,"b_ecal_bar_r_layer_14/F");
    outtree->Branch("ecal_shower_nhit_max_layer",&b_ecal_shower_nhit_max_layer,"b_ecal_shower_nhit_max_layer/F");
    outtree->Branch("ecal_shower_nhit_start_layer",&b_ecal_shower_nhit_start_layer,"b_ecal_shower_nhit_start_layer/F");
    outtree->Branch("ecal_shower_nhit_end_layer",&b_ecal_shower_nhit_end_layer,"b_ecal_shower_nhit_end_layer/F");
    outtree->Branch("ecal_shower_nhit_start_10_layer",&b_ecal_shower_nhit_start_10_layer,"b_ecal_shower_nhit_start_10_layer/F");
    outtree->Branch("ecal_shower_nhit_end_10_layer",&b_ecal_shower_nhit_end_10_layer,"b_ecal_shower_nhit_end_10_layer/F");
    outtree->Branch("ecal_shower_nhit_average",&b_ecal_shower_nhit_average,"b_ecal_shower_nhit_average/F");
    outtree->Branch("ecal_shower_nhit_max",&b_ecal_shower_nhit_max,"b_ecal_shower_nhit_max/F");
    outtree->Branch("ecal_shower_sume_max_layer",&b_ecal_shower_sume_max_layer,"b_ecal_shower_sume_max_layer/F");
    outtree->Branch("ecal_shower_sume_start_layer",&b_ecal_shower_sume_start_layer,"b_ecal_shower_sume_start_layer/F");
    outtree->Branch("ecal_shower_sume_end_layer",&b_ecal_shower_sume_end_layer,"b_ecal_shower_sume_end_layer/F");
    outtree->Branch("ecal_shower_sume_start_10_layer",&b_ecal_shower_sume_start_10_layer,"b_ecal_shower_sume_start_10_layer/F");
    outtree->Branch("ecal_shower_sume_end_10_layer",&b_ecal_shower_sume_end_10_layer,"b_ecal_shower_sume_end_10_layer/F");
    outtree->Branch("ecal_shower_sume_average",&b_ecal_shower_sume_average,"b_ecal_shower_sume_average/F");
    outtree->Branch("ecal_shower_sume_max",&b_ecal_shower_sume_max,"b_ecal_shower_sume_max/F");
    outtree->Branch("ecal_shower_weighte_max_layer",&b_ecal_shower_weighte_max_layer,"b_ecal_shower_weighte_max_layer/F");
    outtree->Branch("ecal_shower_weighte_start_layer",&b_ecal_shower_weighte_start_layer,"b_ecal_shower_weighte_start_layer/F");
    outtree->Branch("ecal_shower_weighte_end_layer",&b_ecal_shower_weighte_end_layer,"b_ecal_shower_weighte_end_layer/F");
    outtree->Branch("ecal_shower_weighte_start_10_layer",&b_ecal_shower_weighte_start_10_layer,"b_ecal_shower_weighte_start_10_layer/F");
    outtree->Branch("ecal_shower_weighte_end_10_layer",&b_ecal_shower_weighte_end_10_layer,"b_ecal_shower_weighte_end_10_layer/F");
    outtree->Branch("ecal_shower_weighte_average",&b_ecal_shower_weighte_average,"b_ecal_shower_weighte_average/F");
    outtree->Branch("ecal_shower_weighte_max",&b_ecal_shower_weighte_max,"b_ecal_shower_weighte_max/F");
    outtree->Branch("ecal_nhit_layer_0",&b_ecal_nhit_layer_0,"b_ecal_nhit_layer_0/F");
    outtree->Branch("ecal_nhit_layer_1",&b_ecal_nhit_layer_1,"b_ecal_nhit_layer_1/F");
    outtree->Branch("ecal_nhit_layer_2",&b_ecal_nhit_layer_2,"b_ecal_nhit_layer_2/F");
    outtree->Branch("ecal_nhit_layer_3",&b_ecal_nhit_layer_3,"b_ecal_nhit_layer_3/F");
    outtree->Branch("ecal_nhit_layer_4",&b_ecal_nhit_layer_4,"b_ecal_nhit_layer_4/F");
    outtree->Branch("ecal_nhit_layer_5",&b_ecal_nhit_layer_5,"b_ecal_nhit_layer_5/F");
    outtree->Branch("ecal_nhit_layer_6",&b_ecal_nhit_layer_6,"b_ecal_nhit_layer_6/F");
    outtree->Branch("ecal_nhit_layer_7",&b_ecal_nhit_layer_7,"b_ecal_nhit_layer_7/F");
    outtree->Branch("ecal_nhit_layer_8",&b_ecal_nhit_layer_8,"b_ecal_nhit_layer_8/F");
    outtree->Branch("ecal_nhit_layer_9",&b_ecal_nhit_layer_9,"b_ecal_nhit_layer_9/F");
    outtree->Branch("ecal_nhit_layer_10",&b_ecal_nhit_layer_10,"b_ecal_nhit_layer_10/F");
    outtree->Branch("ecal_nhit_layer_11",&b_ecal_nhit_layer_11,"b_ecal_nhit_layer_11/F");
    outtree->Branch("ecal_nhit_layer_12",&b_ecal_nhit_layer_12,"b_ecal_nhit_layer_12/F");
    outtree->Branch("ecal_nhit_layer_13",&b_ecal_nhit_layer_13,"b_ecal_nhit_layer_13/F");
    outtree->Branch("ecal_nhit_layer_14",&b_ecal_nhit_layer_14,"b_ecal_nhit_layer_14/F");
    outtree->Branch("ecal_nhit_layer_n_0",&b_ecal_nhit_layer_n_0,"b_ecal_nhit_layer_n_0/F");
    outtree->Branch("ecal_nhit_layer_n_1",&b_ecal_nhit_layer_n_1,"b_ecal_nhit_layer_n_1/F");
    outtree->Branch("ecal_nhit_layer_n_2",&b_ecal_nhit_layer_n_2,"b_ecal_nhit_layer_n_2/F");
    outtree->Branch("ecal_nhit_layer_n_3",&b_ecal_nhit_layer_n_3,"b_ecal_nhit_layer_n_3/F");
    outtree->Branch("ecal_nhit_layer_n_4",&b_ecal_nhit_layer_n_4,"b_ecal_nhit_layer_n_4/F");
    outtree->Branch("ecal_nhit_layer_n_5",&b_ecal_nhit_layer_n_5,"b_ecal_nhit_layer_n_5/F");
    outtree->Branch("ecal_nhit_layer_n_6",&b_ecal_nhit_layer_n_6,"b_ecal_nhit_layer_n_6/F");
    outtree->Branch("ecal_nhit_layer_n_7",&b_ecal_nhit_layer_n_7,"b_ecal_nhit_layer_n_7/F");
    outtree->Branch("ecal_nhit_layer_n_8",&b_ecal_nhit_layer_n_8,"b_ecal_nhit_layer_n_8/F");
    outtree->Branch("ecal_nhit_layer_n_9",&b_ecal_nhit_layer_n_9,"b_ecal_nhit_layer_n_9/F");
    outtree->Branch("ecal_nhit_layer_n_10",&b_ecal_nhit_layer_n_10,"b_ecal_nhit_layer_n_10/F");
    outtree->Branch("ecal_nhit_layer_n_11",&b_ecal_nhit_layer_n_11,"b_ecal_nhit_layer_n_11/F");
    outtree->Branch("ecal_nhit_layer_n_12",&b_ecal_nhit_layer_n_12,"b_ecal_nhit_layer_n_12/F");
    outtree->Branch("ecal_nhit_layer_n_13",&b_ecal_nhit_layer_n_13,"b_ecal_nhit_layer_n_13/F");
    outtree->Branch("ecal_nhit_layer_n_14",&b_ecal_nhit_layer_n_14,"b_ecal_nhit_layer_n_14/F");
    outtree->Branch("ecal_sume_layer_0",&b_ecal_sume_layer_0,"b_ecal_sume_layer_0/F");
    outtree->Branch("ecal_sume_layer_1",&b_ecal_sume_layer_1,"b_ecal_sume_layer_1/F");
    outtree->Branch("ecal_sume_layer_2",&b_ecal_sume_layer_2,"b_ecal_sume_layer_2/F");
    outtree->Branch("ecal_sume_layer_3",&b_ecal_sume_layer_3,"b_ecal_sume_layer_3/F");
    outtree->Branch("ecal_sume_layer_4",&b_ecal_sume_layer_4,"b_ecal_sume_layer_4/F");
    outtree->Branch("ecal_sume_layer_5",&b_ecal_sume_layer_5,"b_ecal_sume_layer_5/F");
    outtree->Branch("ecal_sume_layer_6",&b_ecal_sume_layer_6,"b_ecal_sume_layer_6/F");
    outtree->Branch("ecal_sume_layer_7",&b_ecal_sume_layer_7,"b_ecal_sume_layer_7/F");
    outtree->Branch("ecal_sume_layer_8",&b_ecal_sume_layer_8,"b_ecal_sume_layer_8/F");
    outtree->Branch("ecal_sume_layer_9",&b_ecal_sume_layer_9,"b_ecal_sume_layer_9/F");
    outtree->Branch("ecal_sume_layer_10",&b_ecal_sume_layer_10,"b_ecal_sume_layer_10/F");
    outtree->Branch("ecal_sume_layer_11",&b_ecal_sume_layer_11,"b_ecal_sume_layer_11/F");
    outtree->Branch("ecal_sume_layer_12",&b_ecal_sume_layer_12,"b_ecal_sume_layer_12/F");
    outtree->Branch("ecal_sume_layer_13",&b_ecal_sume_layer_13,"b_ecal_sume_layer_13/F");
    outtree->Branch("ecal_sume_layer_14",&b_ecal_sume_layer_14,"b_ecal_sume_layer_14/F");
    outtree->Branch("ecal_sume_layer_n_0",&b_ecal_sume_layer_n_0,"b_ecal_sume_layer_n_0/F");
    outtree->Branch("ecal_sume_layer_n_1",&b_ecal_sume_layer_n_1,"b_ecal_sume_layer_n_1/F");
    outtree->Branch("ecal_sume_layer_n_2",&b_ecal_sume_layer_n_2,"b_ecal_sume_layer_n_2/F");
    outtree->Branch("ecal_sume_layer_n_3",&b_ecal_sume_layer_n_3,"b_ecal_sume_layer_n_3/F");
    outtree->Branch("ecal_sume_layer_n_4",&b_ecal_sume_layer_n_4,"b_ecal_sume_layer_n_4/F");
    outtree->Branch("ecal_sume_layer_n_5",&b_ecal_sume_layer_n_5,"b_ecal_sume_layer_n_5/F");
    outtree->Branch("ecal_sume_layer_n_6",&b_ecal_sume_layer_n_6,"b_ecal_sume_layer_n_6/F");
    outtree->Branch("ecal_sume_layer_n_7",&b_ecal_sume_layer_n_7,"b_ecal_sume_layer_n_7/F");
    outtree->Branch("ecal_sume_layer_n_8",&b_ecal_sume_layer_n_8,"b_ecal_sume_layer_n_8/F");
    outtree->Branch("ecal_sume_layer_n_9",&b_ecal_sume_layer_n_9,"b_ecal_sume_layer_n_9/F");
    outtree->Branch("ecal_sume_layer_n_10",&b_ecal_sume_layer_n_10,"b_ecal_sume_layer_n_10/F");
    outtree->Branch("ecal_sume_layer_n_11",&b_ecal_sume_layer_n_11,"b_ecal_sume_layer_n_11/F");
    outtree->Branch("ecal_sume_layer_n_12",&b_ecal_sume_layer_n_12,"b_ecal_sume_layer_n_12/F");
    outtree->Branch("ecal_sume_layer_n_13",&b_ecal_sume_layer_n_13,"b_ecal_sume_layer_n_13/F");
    outtree->Branch("ecal_sume_layer_n_14",&b_ecal_sume_layer_n_14,"b_ecal_sume_layer_n_14/F");
    outtree->Branch("ecal_weighte_layer_0",&b_ecal_weighte_layer_0,"b_ecal_weighte_layer_0/F");
    outtree->Branch("ecal_weighte_layer_1",&b_ecal_weighte_layer_1,"b_ecal_weighte_layer_1/F");
    outtree->Branch("ecal_weighte_layer_2",&b_ecal_weighte_layer_2,"b_ecal_weighte_layer_2/F");
    outtree->Branch("ecal_weighte_layer_3",&b_ecal_weighte_layer_3,"b_ecal_weighte_layer_3/F");
    outtree->Branch("ecal_weighte_layer_4",&b_ecal_weighte_layer_4,"b_ecal_weighte_layer_4/F");
    outtree->Branch("ecal_weighte_layer_5",&b_ecal_weighte_layer_5,"b_ecal_weighte_layer_5/F");
    outtree->Branch("ecal_weighte_layer_6",&b_ecal_weighte_layer_6,"b_ecal_weighte_layer_6/F");
    outtree->Branch("ecal_weighte_layer_7",&b_ecal_weighte_layer_7,"b_ecal_weighte_layer_7/F");
    outtree->Branch("ecal_weighte_layer_8",&b_ecal_weighte_layer_8,"b_ecal_weighte_layer_8/F");
    outtree->Branch("ecal_weighte_layer_9",&b_ecal_weighte_layer_9,"b_ecal_weighte_layer_9/F");
    outtree->Branch("ecal_weighte_layer_10",&b_ecal_weighte_layer_10,"b_ecal_weighte_layer_10/F");
    outtree->Branch("ecal_weighte_layer_11",&b_ecal_weighte_layer_11,"b_ecal_weighte_layer_11/F");
    outtree->Branch("ecal_weighte_layer_12",&b_ecal_weighte_layer_12,"b_ecal_weighte_layer_12/F");
    outtree->Branch("ecal_weighte_layer_13",&b_ecal_weighte_layer_13,"b_ecal_weighte_layer_13/F");
    outtree->Branch("ecal_weighte_layer_14",&b_ecal_weighte_layer_14,"b_ecal_weighte_layer_14/F");
    outtree->Branch("ecal_weighte_layer_n_0",&b_ecal_weighte_layer_n_0,"b_ecal_weighte_layer_n_0/F");
    outtree->Branch("ecal_weighte_layer_n_1",&b_ecal_weighte_layer_n_1,"b_ecal_weighte_layer_n_1/F");
    outtree->Branch("ecal_weighte_layer_n_2",&b_ecal_weighte_layer_n_2,"b_ecal_weighte_layer_n_2/F");
    outtree->Branch("ecal_weighte_layer_n_3",&b_ecal_weighte_layer_n_3,"b_ecal_weighte_layer_n_3/F");
    outtree->Branch("ecal_weighte_layer_n_4",&b_ecal_weighte_layer_n_4,"b_ecal_weighte_layer_n_4/F");
    outtree->Branch("ecal_weighte_layer_n_5",&b_ecal_weighte_layer_n_5,"b_ecal_weighte_layer_n_5/F");
    outtree->Branch("ecal_weighte_layer_n_6",&b_ecal_weighte_layer_n_6,"b_ecal_weighte_layer_n_6/F");
    outtree->Branch("ecal_weighte_layer_n_7",&b_ecal_weighte_layer_n_7,"b_ecal_weighte_layer_n_7/F");
    outtree->Branch("ecal_weighte_layer_n_8",&b_ecal_weighte_layer_n_8,"b_ecal_weighte_layer_n_8/F");
    outtree->Branch("ecal_weighte_layer_n_9",&b_ecal_weighte_layer_n_9,"b_ecal_weighte_layer_n_9/F");
    outtree->Branch("ecal_weighte_layer_n_10",&b_ecal_weighte_layer_n_10,"b_ecal_weighte_layer_n_10/F");
    outtree->Branch("ecal_weighte_layer_n_11",&b_ecal_weighte_layer_n_11,"b_ecal_weighte_layer_n_11/F");
    outtree->Branch("ecal_weighte_layer_n_12",&b_ecal_weighte_layer_n_12,"b_ecal_weighte_layer_n_12/F");
    outtree->Branch("ecal_weighte_layer_n_13",&b_ecal_weighte_layer_n_13,"b_ecal_weighte_layer_n_13/F");
    outtree->Branch("ecal_weighte_layer_n_14",&b_ecal_weighte_layer_n_14,"b_ecal_weighte_layer_n_14/F");
    outtree->Branch("ecal_sume_layer_0",&b_ecal_sume_layer_0,"b_ecal_sume_layer_0/F");
    outtree->Branch("ecal_sume_layer_1",&b_ecal_sume_layer_1,"b_ecal_sume_layer_1/F");
    outtree->Branch("ecal_sume_layer_2",&b_ecal_sume_layer_2,"b_ecal_sume_layer_2/F");
    outtree->Branch("ecal_sume_layer_3",&b_ecal_sume_layer_3,"b_ecal_sume_layer_3/F");
    outtree->Branch("ecal_sume_layer_4",&b_ecal_sume_layer_4,"b_ecal_sume_layer_4/F");
    outtree->Branch("ecal_sume_layer_5",&b_ecal_sume_layer_5,"b_ecal_sume_layer_5/F");
    outtree->Branch("ecal_sume_layer_6",&b_ecal_sume_layer_6,"b_ecal_sume_layer_6/F");
    outtree->Branch("ecal_sume_layer_7",&b_ecal_sume_layer_7,"b_ecal_sume_layer_7/F");
    outtree->Branch("ecal_sume_layer_8",&b_ecal_sume_layer_8,"b_ecal_sume_layer_8/F");
    outtree->Branch("ecal_sume_layer_9",&b_ecal_sume_layer_9,"b_ecal_sume_layer_9/F");
    outtree->Branch("ecal_sume_layer_10",&b_ecal_sume_layer_10,"b_ecal_sume_layer_10/F");
    outtree->Branch("ecal_sume_layer_11",&b_ecal_sume_layer_11,"b_ecal_sume_layer_11/F");
    outtree->Branch("ecal_sume_layer_12",&b_ecal_sume_layer_12,"b_ecal_sume_layer_12/F");
    outtree->Branch("ecal_sume_layer_13",&b_ecal_sume_layer_13,"b_ecal_sume_layer_13/F");
    outtree->Branch("ecal_sume_layer_14",&b_ecal_sume_layer_14,"b_ecal_sume_layer_14/F");
    outtree->Branch("ecal_sume_layer_n_0",&b_ecal_sume_layer_n_0,"b_ecal_sume_layer_n_0/F");
    outtree->Branch("ecal_sume_layer_n_1",&b_ecal_sume_layer_n_1,"b_ecal_sume_layer_n_1/F");
    outtree->Branch("ecal_sume_layer_n_2",&b_ecal_sume_layer_n_2,"b_ecal_sume_layer_n_2/F");
    outtree->Branch("ecal_sume_layer_n_3",&b_ecal_sume_layer_n_3,"b_ecal_sume_layer_n_3/F");
    outtree->Branch("ecal_sume_layer_n_4",&b_ecal_sume_layer_n_4,"b_ecal_sume_layer_n_4/F");
    outtree->Branch("ecal_sume_layer_n_5",&b_ecal_sume_layer_n_5,"b_ecal_sume_layer_n_5/F");
    outtree->Branch("ecal_sume_layer_n_6",&b_ecal_sume_layer_n_6,"b_ecal_sume_layer_n_6/F");
    outtree->Branch("ecal_sume_layer_n_7",&b_ecal_sume_layer_n_7,"b_ecal_sume_layer_n_7/F");
    outtree->Branch("ecal_sume_layer_n_8",&b_ecal_sume_layer_n_8,"b_ecal_sume_layer_n_8/F");
    outtree->Branch("ecal_sume_layer_n_9",&b_ecal_sume_layer_n_9,"b_ecal_sume_layer_n_9/F");
    outtree->Branch("ecal_sume_layer_n_10",&b_ecal_sume_layer_n_10,"b_ecal_sume_layer_n_10/F");
    outtree->Branch("ecal_sume_layer_n_11",&b_ecal_sume_layer_n_11,"b_ecal_sume_layer_n_11/F");
    outtree->Branch("ecal_sume_layer_n_12",&b_ecal_sume_layer_n_12,"b_ecal_sume_layer_n_12/F");
    outtree->Branch("ecal_sume_layer_n_13",&b_ecal_sume_layer_n_13,"b_ecal_sume_layer_n_13/F");
    outtree->Branch("ecal_sume_layer_n_14",&b_ecal_sume_layer_n_14,"b_ecal_sume_layer_n_14/F");

    outtree->Branch("hcal_interaction",&b_hcal_interaction,"b_hcal_interaction/F");
    outtree->Branch("hcal_nhit",&b_hcal_nhit,"b_hcal_nhit/F");
    outtree->Branch("hcal_sume",&b_hcal_sume,"b_hcal_sume/F");
    outtree->Branch("hcal_weighte",&b_hcal_weighte,"b_hcal_weighte/F");
    outtree->Branch("hcal_mol",&b_hcal_mol,"b_hcal_mol/F");
    outtree->Branch("hcal_MIP_Likeness",&b_hcal_MIP_Likeness,"b_hcal_MIP_Likeness/F");
    outtree->Branch("hcal_hits_max_distance",&b_hcal_hits_max_distance,"b_hcal_hits_max_distance/F");
    outtree->Branch("hcal_radius90_layer_0",&b_hcal_radius90_layer_0,"b_hcal_radius90_layer_0/F");
    outtree->Branch("hcal_radius90_layer_1",&b_hcal_radius90_layer_1,"b_hcal_radius90_layer_1/F");
    outtree->Branch("hcal_radius90_layer_2",&b_hcal_radius90_layer_2,"b_hcal_radius90_layer_2/F");
    outtree->Branch("hcal_radius90_layer_3",&b_hcal_radius90_layer_3,"b_hcal_radius90_layer_3/F");
    outtree->Branch("hcal_radius90_layer_4",&b_hcal_radius90_layer_4,"b_hcal_radius90_layer_4/F");
    outtree->Branch("hcal_radius90_layer_5",&b_hcal_radius90_layer_5,"b_hcal_radius90_layer_5/F");
    outtree->Branch("hcal_radius90_layer_6",&b_hcal_radius90_layer_6,"b_hcal_radius90_layer_6/F");
    outtree->Branch("hcal_radius90_layer_7",&b_hcal_radius90_layer_7,"b_hcal_radius90_layer_7/F");
    outtree->Branch("hcal_radius90_layer_8",&b_hcal_radius90_layer_8,"b_hcal_radius90_layer_8/F");
    outtree->Branch("hcal_radius90_layer_9",&b_hcal_radius90_layer_9,"b_hcal_radius90_layer_9/F");
    outtree->Branch("hcal_radius90_layer_10",&b_hcal_radius90_layer_10,"b_hcal_radius90_layer_10/F");
    outtree->Branch("hcal_radius90_layer_11",&b_hcal_radius90_layer_11,"b_hcal_radius90_layer_11/F");
    outtree->Branch("hcal_radius90_layer_12",&b_hcal_radius90_layer_12,"b_hcal_radius90_layer_12/F");
    outtree->Branch("hcal_radius90_layer_13",&b_hcal_radius90_layer_13,"b_hcal_radius90_layer_13/F");
    outtree->Branch("hcal_radius90_layer_14",&b_hcal_radius90_layer_14,"b_hcal_radius90_layer_14/F");
    outtree->Branch("hcal_bar_x",&b_hcal_bar_x,"b_hcal_bar_x/F");
    outtree->Branch("hcal_bar_y",&b_hcal_bar_y,"b_hcal_bar_y/F");
    outtree->Branch("hcal_bar_z",&b_hcal_bar_z,"b_hcal_bar_z/F");
    outtree->Branch("hcal_bar_r",&b_hcal_bar_r,"b_hcal_bar_r/F");
    outtree->Branch("hcal_bar_x_layer_0",&b_hcal_bar_x_layer_0,"b_hcal_bar_x_layer_0/F");
    outtree->Branch("hcal_bar_x_layer_1",&b_hcal_bar_x_layer_1,"b_hcal_bar_x_layer_1/F");
    outtree->Branch("hcal_bar_x_layer_2",&b_hcal_bar_x_layer_2,"b_hcal_bar_x_layer_2/F");
    outtree->Branch("hcal_bar_x_layer_3",&b_hcal_bar_x_layer_3,"b_hcal_bar_x_layer_3/F");
    outtree->Branch("hcal_bar_x_layer_4",&b_hcal_bar_x_layer_4,"b_hcal_bar_x_layer_4/F");
    outtree->Branch("hcal_bar_x_layer_5",&b_hcal_bar_x_layer_5,"b_hcal_bar_x_layer_5/F");
    outtree->Branch("hcal_bar_x_layer_6",&b_hcal_bar_x_layer_6,"b_hcal_bar_x_layer_6/F");
    outtree->Branch("hcal_bar_x_layer_7",&b_hcal_bar_x_layer_7,"b_hcal_bar_x_layer_7/F");
    outtree->Branch("hcal_bar_x_layer_8",&b_hcal_bar_x_layer_8,"b_hcal_bar_x_layer_8/F");
    outtree->Branch("hcal_bar_x_layer_9",&b_hcal_bar_x_layer_9,"b_hcal_bar_x_layer_9/F");
    outtree->Branch("hcal_bar_x_layer_10",&b_hcal_bar_x_layer_10,"b_hcal_bar_x_layer_10/F");
    outtree->Branch("hcal_bar_x_layer_11",&b_hcal_bar_x_layer_11,"b_hcal_bar_x_layer_11/F");
    outtree->Branch("hcal_bar_x_layer_12",&b_hcal_bar_x_layer_12,"b_hcal_bar_x_layer_12/F");
    outtree->Branch("hcal_bar_x_layer_13",&b_hcal_bar_x_layer_13,"b_hcal_bar_x_layer_13/F");
    outtree->Branch("hcal_bar_x_layer_14",&b_hcal_bar_x_layer_14,"b_hcal_bar_x_layer_14/F");
    outtree->Branch("hcal_bar_y_layer_0",&b_hcal_bar_y_layer_0,"b_hcal_bar_y_layer_0/F");
    outtree->Branch("hcal_bar_y_layer_1",&b_hcal_bar_y_layer_1,"b_hcal_bar_y_layer_1/F");
    outtree->Branch("hcal_bar_y_layer_2",&b_hcal_bar_y_layer_2,"b_hcal_bar_y_layer_2/F");
    outtree->Branch("hcal_bar_y_layer_3",&b_hcal_bar_y_layer_3,"b_hcal_bar_y_layer_3/F");
    outtree->Branch("hcal_bar_y_layer_4",&b_hcal_bar_y_layer_4,"b_hcal_bar_y_layer_4/F");
    outtree->Branch("hcal_bar_y_layer_5",&b_hcal_bar_y_layer_5,"b_hcal_bar_y_layer_5/F");
    outtree->Branch("hcal_bar_y_layer_6",&b_hcal_bar_y_layer_6,"b_hcal_bar_y_layer_6/F");
    outtree->Branch("hcal_bar_y_layer_7",&b_hcal_bar_y_layer_7,"b_hcal_bar_y_layer_7/F");
    outtree->Branch("hcal_bar_y_layer_8",&b_hcal_bar_y_layer_8,"b_hcal_bar_y_layer_8/F");
    outtree->Branch("hcal_bar_y_layer_9",&b_hcal_bar_y_layer_9,"b_hcal_bar_y_layer_9/F");
    outtree->Branch("hcal_bar_y_layer_10",&b_hcal_bar_y_layer_10,"b_hcal_bar_y_layer_10/F");
    outtree->Branch("hcal_bar_y_layer_11",&b_hcal_bar_y_layer_11,"b_hcal_bar_y_layer_11/F");
    outtree->Branch("hcal_bar_y_layer_12",&b_hcal_bar_y_layer_12,"b_hcal_bar_y_layer_12/F");
    outtree->Branch("hcal_bar_y_layer_13",&b_hcal_bar_y_layer_13,"b_hcal_bar_y_layer_13/F");
    outtree->Branch("hcal_bar_y_layer_14",&b_hcal_bar_y_layer_14,"b_hcal_bar_y_layer_14/F");
    outtree->Branch("hcal_bar_r_layer_0",&b_hcal_bar_r_layer_0,"b_hcal_bar_r_layer_0/F");
    outtree->Branch("hcal_bar_r_layer_1",&b_hcal_bar_r_layer_1,"b_hcal_bar_r_layer_1/F");
    outtree->Branch("hcal_bar_r_layer_2",&b_hcal_bar_r_layer_2,"b_hcal_bar_r_layer_2/F");
    outtree->Branch("hcal_bar_r_layer_3",&b_hcal_bar_r_layer_3,"b_hcal_bar_r_layer_3/F");
    outtree->Branch("hcal_bar_r_layer_4",&b_hcal_bar_r_layer_4,"b_hcal_bar_r_layer_4/F");
    outtree->Branch("hcal_bar_r_layer_5",&b_hcal_bar_r_layer_5,"b_hcal_bar_r_layer_5/F");
    outtree->Branch("hcal_bar_r_layer_6",&b_hcal_bar_r_layer_6,"b_hcal_bar_r_layer_6/F");
    outtree->Branch("hcal_bar_r_layer_7",&b_hcal_bar_r_layer_7,"b_hcal_bar_r_layer_7/F");
    outtree->Branch("hcal_bar_r_layer_8",&b_hcal_bar_r_layer_8,"b_hcal_bar_r_layer_8/F");
    outtree->Branch("hcal_bar_r_layer_9",&b_hcal_bar_r_layer_9,"b_hcal_bar_r_layer_9/F");
    outtree->Branch("hcal_bar_r_layer_10",&b_hcal_bar_r_layer_10,"b_hcal_bar_r_layer_10/F");
    outtree->Branch("hcal_bar_r_layer_11",&b_hcal_bar_r_layer_11,"b_hcal_bar_r_layer_11/F");
    outtree->Branch("hcal_bar_r_layer_12",&b_hcal_bar_r_layer_12,"b_hcal_bar_r_layer_12/F");
    outtree->Branch("hcal_bar_r_layer_13",&b_hcal_bar_r_layer_13,"b_hcal_bar_r_layer_13/F");
    outtree->Branch("hcal_bar_r_layer_14",&b_hcal_bar_r_layer_14,"b_hcal_bar_r_layer_14/F");
    outtree->Branch("hcal_shower_nhit_max_layer",&b_hcal_shower_nhit_max_layer,"b_hcal_shower_nhit_max_layer/F");
    outtree->Branch("hcal_shower_nhit_start_layer",&b_hcal_shower_nhit_start_layer,"b_hcal_shower_nhit_start_layer/F");
    outtree->Branch("hcal_shower_nhit_end_layer",&b_hcal_shower_nhit_end_layer,"b_hcal_shower_nhit_end_layer/F");
    outtree->Branch("hcal_shower_nhit_start_10_layer",&b_hcal_shower_nhit_start_10_layer,"b_hcal_shower_nhit_start_10_layer/F");
    outtree->Branch("hcal_shower_nhit_end_10_layer",&b_hcal_shower_nhit_end_10_layer,"b_hcal_shower_nhit_end_10_layer/F");
    outtree->Branch("hcal_shower_nhit_average",&b_hcal_shower_nhit_average,"b_hcal_shower_nhit_average/F");
    outtree->Branch("hcal_shower_nhit_max",&b_hcal_shower_nhit_max,"b_hcal_shower_nhit_max/F");
    outtree->Branch("hcal_shower_sume_max_layer",&b_hcal_shower_sume_max_layer,"b_hcal_shower_sume_max_layer/F");
    outtree->Branch("hcal_shower_sume_start_layer",&b_hcal_shower_sume_start_layer,"b_hcal_shower_sume_start_layer/F");
    outtree->Branch("hcal_shower_sume_end_layer",&b_hcal_shower_sume_end_layer,"b_hcal_shower_sume_end_layer/F");
    outtree->Branch("hcal_shower_sume_start_10_layer",&b_hcal_shower_sume_start_10_layer,"b_hcal_shower_sume_start_10_layer/F");
    outtree->Branch("hcal_shower_sume_end_10_layer",&b_hcal_shower_sume_end_10_layer,"b_hcal_shower_sume_end_10_layer/F");
    outtree->Branch("hcal_shower_sume_average",&b_hcal_shower_sume_average,"b_hcal_shower_sume_average/F");
    outtree->Branch("hcal_shower_sume_max",&b_hcal_shower_sume_max,"b_hcal_shower_sume_max/F");
    outtree->Branch("hcal_shower_weighte_max_layer",&b_hcal_shower_weighte_max_layer,"b_hcal_shower_weighte_max_layer/F");
    outtree->Branch("hcal_shower_weighte_start_layer",&b_hcal_shower_weighte_start_layer,"b_hcal_shower_weighte_start_layer/F");
    outtree->Branch("hcal_shower_weighte_end_layer",&b_hcal_shower_weighte_end_layer,"b_hcal_shower_weighte_end_layer/F");
    outtree->Branch("hcal_shower_weighte_start_10_layer",&b_hcal_shower_weighte_start_10_layer,"b_hcal_shower_weighte_start_10_layer/F");
    outtree->Branch("hcal_shower_weighte_end_10_layer",&b_hcal_shower_weighte_end_10_layer,"b_hcal_shower_weighte_end_10_layer/F");
    outtree->Branch("hcal_shower_weighte_average",&b_hcal_shower_weighte_average,"b_hcal_shower_weighte_average/F");
    outtree->Branch("hcal_shower_weighte_max",&b_hcal_shower_weighte_max,"b_hcal_shower_weighte_max/F");
    outtree->Branch("hcal_nhit_layer_0", &b_hcal_nhit_layer_0, "b_hcal_nhit_layer_0/F");
    outtree->Branch("hcal_nhit_layer_1", &b_hcal_nhit_layer_1, "b_hcal_nhit_layer_1/F");
    outtree->Branch("hcal_nhit_layer_2", &b_hcal_nhit_layer_2, "b_hcal_nhit_layer_2/F");
    outtree->Branch("hcal_nhit_layer_3", &b_hcal_nhit_layer_3, "b_hcal_nhit_layer_3/F");
    outtree->Branch("hcal_nhit_layer_4", &b_hcal_nhit_layer_4, "b_hcal_nhit_layer_4/F");
    outtree->Branch("hcal_nhit_layer_5", &b_hcal_nhit_layer_5, "b_hcal_nhit_layer_5/F");
    outtree->Branch("hcal_nhit_layer_6", &b_hcal_nhit_layer_6, "b_hcal_nhit_layer_6/F");
    outtree->Branch("hcal_nhit_layer_7", &b_hcal_nhit_layer_7, "b_hcal_nhit_layer_7/F");
    outtree->Branch("hcal_nhit_layer_8", &b_hcal_nhit_layer_8, "b_hcal_nhit_layer_8/F");
    outtree->Branch("hcal_nhit_layer_9", &b_hcal_nhit_layer_9, "b_hcal_nhit_layer_9/F");
    outtree->Branch("hcal_nhit_layer_10", &b_hcal_nhit_layer_10, "b_hcal_nhit_layer_10/F");
    outtree->Branch("hcal_nhit_layer_11", &b_hcal_nhit_layer_11, "b_hcal_nhit_layer_11/F");
    outtree->Branch("hcal_nhit_layer_12", &b_hcal_nhit_layer_12, "b_hcal_nhit_layer_12/F");
    outtree->Branch("hcal_nhit_layer_13", &b_hcal_nhit_layer_13, "b_hcal_nhit_layer_13/F");
    outtree->Branch("hcal_nhit_layer_14", &b_hcal_nhit_layer_14, "b_hcal_nhit_layer_14/F");
    outtree->Branch("hcal_nhit_layer_15", &b_hcal_nhit_layer_15, "b_hcal_nhit_layer_15/F");
    outtree->Branch("hcal_nhit_layer_16", &b_hcal_nhit_layer_16, "b_hcal_nhit_layer_16/F");
    outtree->Branch("hcal_nhit_layer_17", &b_hcal_nhit_layer_17, "b_hcal_nhit_layer_17/F");
    outtree->Branch("hcal_nhit_layer_18", &b_hcal_nhit_layer_18, "b_hcal_nhit_layer_18/F");
    outtree->Branch("hcal_nhit_layer_19", &b_hcal_nhit_layer_19, "b_hcal_nhit_layer_19/F");
    outtree->Branch("hcal_nhit_layer_20", &b_hcal_nhit_layer_20, "b_hcal_nhit_layer_20/F");
    outtree->Branch("hcal_nhit_layer_21", &b_hcal_nhit_layer_21, "b_hcal_nhit_layer_21/F");
    outtree->Branch("hcal_nhit_layer_22", &b_hcal_nhit_layer_22, "b_hcal_nhit_layer_22/F");
    outtree->Branch("hcal_nhit_layer_23", &b_hcal_nhit_layer_23, "b_hcal_nhit_layer_23/F");
    outtree->Branch("hcal_nhit_layer_24", &b_hcal_nhit_layer_24, "b_hcal_nhit_layer_24/F");
    outtree->Branch("hcal_nhit_layer_25", &b_hcal_nhit_layer_25, "b_hcal_nhit_layer_25/F");
    outtree->Branch("hcal_nhit_layer_26", &b_hcal_nhit_layer_26, "b_hcal_nhit_layer_26/F");
    outtree->Branch("hcal_nhit_layer_27", &b_hcal_nhit_layer_27, "b_hcal_nhit_layer_27/F");
    outtree->Branch("hcal_nhit_layer_28", &b_hcal_nhit_layer_28, "b_hcal_nhit_layer_28/F");
    outtree->Branch("hcal_nhit_layer_29", &b_hcal_nhit_layer_29, "b_hcal_nhit_layer_29/F");
    outtree->Branch("hcal_nhit_layer_30", &b_hcal_nhit_layer_30, "b_hcal_nhit_layer_30/F");
    outtree->Branch("hcal_nhit_layer_31", &b_hcal_nhit_layer_31, "b_hcal_nhit_layer_31/F");
    outtree->Branch("hcal_nhit_layer_32", &b_hcal_nhit_layer_32, "b_hcal_nhit_layer_32/F");
    outtree->Branch("hcal_nhit_layer_33", &b_hcal_nhit_layer_33, "b_hcal_nhit_layer_33/F");
    outtree->Branch("hcal_nhit_layer_34", &b_hcal_nhit_layer_34, "b_hcal_nhit_layer_34/F");
    outtree->Branch("hcal_nhit_layer_35", &b_hcal_nhit_layer_35, "b_hcal_nhit_layer_35/F");
    outtree->Branch("hcal_nhit_layer_36", &b_hcal_nhit_layer_36, "b_hcal_nhit_layer_36/F");
    outtree->Branch("hcal_nhit_layer_37", &b_hcal_nhit_layer_37, "b_hcal_nhit_layer_37/F");
    outtree->Branch("hcal_nhit_layer_38", &b_hcal_nhit_layer_38, "b_hcal_nhit_layer_38/F");
    outtree->Branch("hcal_nhit_layer_39", &b_hcal_nhit_layer_39, "b_hcal_nhit_layer_39/F");
    outtree->Branch("hcal_nhit_layer_40", &b_hcal_nhit_layer_40, "b_hcal_nhit_layer_40/F");
    outtree->Branch("hcal_nhit_layer_n_0", &b_hcal_nhit_layer_n_0, "b_hcal_nhit_layer_n_0/F");
    outtree->Branch("hcal_nhit_layer_n_1", &b_hcal_nhit_layer_n_1, "b_hcal_nhit_layer_n_1/F");
    outtree->Branch("hcal_nhit_layer_n_2", &b_hcal_nhit_layer_n_2, "b_hcal_nhit_layer_n_2/F");
    outtree->Branch("hcal_nhit_layer_n_3", &b_hcal_nhit_layer_n_3, "b_hcal_nhit_layer_n_3/F");
    outtree->Branch("hcal_nhit_layer_n_4", &b_hcal_nhit_layer_n_4, "b_hcal_nhit_layer_n_4/F");
    outtree->Branch("hcal_nhit_layer_n_5", &b_hcal_nhit_layer_n_5, "b_hcal_nhit_layer_n_5/F");
    outtree->Branch("hcal_nhit_layer_n_6", &b_hcal_nhit_layer_n_6, "b_hcal_nhit_layer_n_6/F");
    outtree->Branch("hcal_nhit_layer_n_7", &b_hcal_nhit_layer_n_7, "b_hcal_nhit_layer_n_7/F");
    outtree->Branch("hcal_nhit_layer_n_8", &b_hcal_nhit_layer_n_8, "b_hcal_nhit_layer_n_8/F");
    outtree->Branch("hcal_nhit_layer_n_9", &b_hcal_nhit_layer_n_9, "b_hcal_nhit_layer_n_9/F");
    outtree->Branch("hcal_nhit_layer_n_10", &b_hcal_nhit_layer_n_10, "b_hcal_nhit_layer_n_10/F");
    outtree->Branch("hcal_nhit_layer_n_11", &b_hcal_nhit_layer_n_11, "b_hcal_nhit_layer_n_11/F");
    outtree->Branch("hcal_nhit_layer_n_12", &b_hcal_nhit_layer_n_12, "b_hcal_nhit_layer_n_12/F");
    outtree->Branch("hcal_nhit_layer_n_13", &b_hcal_nhit_layer_n_13, "b_hcal_nhit_layer_n_13/F");
    outtree->Branch("hcal_nhit_layer_n_14", &b_hcal_nhit_layer_n_14, "b_hcal_nhit_layer_n_14/F");
    outtree->Branch("hcal_nhit_layer_n_15", &b_hcal_nhit_layer_n_15, "b_hcal_nhit_layer_n_15/F");
    outtree->Branch("hcal_nhit_layer_n_16", &b_hcal_nhit_layer_n_16, "b_hcal_nhit_layer_n_16/F");
    outtree->Branch("hcal_nhit_layer_n_17", &b_hcal_nhit_layer_n_17, "b_hcal_nhit_layer_n_17/F");
    outtree->Branch("hcal_nhit_layer_n_18", &b_hcal_nhit_layer_n_18, "b_hcal_nhit_layer_n_18/F");
    outtree->Branch("hcal_nhit_layer_n_19", &b_hcal_nhit_layer_n_19, "b_hcal_nhit_layer_n_19/F");
    outtree->Branch("hcal_nhit_layer_n_20", &b_hcal_nhit_layer_n_20, "b_hcal_nhit_layer_n_20/F");
    outtree->Branch("hcal_nhit_layer_n_21", &b_hcal_nhit_layer_n_21, "b_hcal_nhit_layer_n_21/F");
    outtree->Branch("hcal_nhit_layer_n_22", &b_hcal_nhit_layer_n_22, "b_hcal_nhit_layer_n_22/F");
    outtree->Branch("hcal_nhit_layer_n_23", &b_hcal_nhit_layer_n_23, "b_hcal_nhit_layer_n_23/F");
    outtree->Branch("hcal_nhit_layer_n_24", &b_hcal_nhit_layer_n_24, "b_hcal_nhit_layer_n_24/F");
    outtree->Branch("hcal_nhit_layer_n_25", &b_hcal_nhit_layer_n_25, "b_hcal_nhit_layer_n_25/F");
    outtree->Branch("hcal_nhit_layer_n_26", &b_hcal_nhit_layer_n_26, "b_hcal_nhit_layer_n_26/F");
    outtree->Branch("hcal_nhit_layer_n_27", &b_hcal_nhit_layer_n_27, "b_hcal_nhit_layer_n_27/F");
    outtree->Branch("hcal_nhit_layer_n_28", &b_hcal_nhit_layer_n_28, "b_hcal_nhit_layer_n_28/F");
    outtree->Branch("hcal_nhit_layer_n_29", &b_hcal_nhit_layer_n_29, "b_hcal_nhit_layer_n_29/F");
    outtree->Branch("hcal_nhit_layer_n_30", &b_hcal_nhit_layer_n_30, "b_hcal_nhit_layer_n_30/F");
    outtree->Branch("hcal_nhit_layer_n_31", &b_hcal_nhit_layer_n_31, "b_hcal_nhit_layer_n_31/F");
    outtree->Branch("hcal_nhit_layer_n_32", &b_hcal_nhit_layer_n_32, "b_hcal_nhit_layer_n_32/F");
    outtree->Branch("hcal_nhit_layer_n_33", &b_hcal_nhit_layer_n_33, "b_hcal_nhit_layer_n_33/F");
    outtree->Branch("hcal_nhit_layer_n_34", &b_hcal_nhit_layer_n_34, "b_hcal_nhit_layer_n_34/F");
    outtree->Branch("hcal_nhit_layer_n_35", &b_hcal_nhit_layer_n_35, "b_hcal_nhit_layer_n_35/F");
    outtree->Branch("hcal_nhit_layer_n_36", &b_hcal_nhit_layer_n_36, "b_hcal_nhit_layer_n_36/F");
    outtree->Branch("hcal_nhit_layer_n_37", &b_hcal_nhit_layer_n_37, "b_hcal_nhit_layer_n_37/F");
    outtree->Branch("hcal_nhit_layer_n_38", &b_hcal_nhit_layer_n_38, "b_hcal_nhit_layer_n_38/F");
    outtree->Branch("hcal_nhit_layer_n_39", &b_hcal_nhit_layer_n_39, "b_hcal_nhit_layer_n_39/F");
    outtree->Branch("hcal_nhit_layer_n_40", &b_hcal_nhit_layer_n_40, "b_hcal_nhit_layer_n_40/F");
    outtree->Branch("hcal_sume_layer_0", &b_hcal_sume_layer_0, "b_hcal_sume_layer_0/F");
    outtree->Branch("hcal_sume_layer_1", &b_hcal_sume_layer_1, "b_hcal_sume_layer_1/F");
    outtree->Branch("hcal_sume_layer_2", &b_hcal_sume_layer_2, "b_hcal_sume_layer_2/F");
    outtree->Branch("hcal_sume_layer_3", &b_hcal_sume_layer_3, "b_hcal_sume_layer_3/F");
    outtree->Branch("hcal_sume_layer_4", &b_hcal_sume_layer_4, "b_hcal_sume_layer_4/F");
    outtree->Branch("hcal_sume_layer_5", &b_hcal_sume_layer_5, "b_hcal_sume_layer_5/F");
    outtree->Branch("hcal_sume_layer_6", &b_hcal_sume_layer_6, "b_hcal_sume_layer_6/F");
    outtree->Branch("hcal_sume_layer_7", &b_hcal_sume_layer_7, "b_hcal_sume_layer_7/F");
    outtree->Branch("hcal_sume_layer_8", &b_hcal_sume_layer_8, "b_hcal_sume_layer_8/F");
    outtree->Branch("hcal_sume_layer_9", &b_hcal_sume_layer_9, "b_hcal_sume_layer_9/F");
    outtree->Branch("hcal_sume_layer_10", &b_hcal_sume_layer_10, "b_hcal_sume_layer_10/F");
    outtree->Branch("hcal_sume_layer_11", &b_hcal_sume_layer_11, "b_hcal_sume_layer_11/F");
    outtree->Branch("hcal_sume_layer_12", &b_hcal_sume_layer_12, "b_hcal_sume_layer_12/F");
    outtree->Branch("hcal_sume_layer_13", &b_hcal_sume_layer_13, "b_hcal_sume_layer_13/F");
    outtree->Branch("hcal_sume_layer_14", &b_hcal_sume_layer_14, "b_hcal_sume_layer_14/F");
    outtree->Branch("hcal_sume_layer_15", &b_hcal_sume_layer_15, "b_hcal_sume_layer_15/F");
    outtree->Branch("hcal_sume_layer_16", &b_hcal_sume_layer_16, "b_hcal_sume_layer_16/F");
    outtree->Branch("hcal_sume_layer_17", &b_hcal_sume_layer_17, "b_hcal_sume_layer_17/F");
    outtree->Branch("hcal_sume_layer_18", &b_hcal_sume_layer_18, "b_hcal_sume_layer_18/F");
    outtree->Branch("hcal_sume_layer_19", &b_hcal_sume_layer_19, "b_hcal_sume_layer_19/F");
    outtree->Branch("hcal_sume_layer_20", &b_hcal_sume_layer_20, "b_hcal_sume_layer_20/F");
    outtree->Branch("hcal_sume_layer_21", &b_hcal_sume_layer_21, "b_hcal_sume_layer_21/F");
    outtree->Branch("hcal_sume_layer_22", &b_hcal_sume_layer_22, "b_hcal_sume_layer_22/F");
    outtree->Branch("hcal_sume_layer_23", &b_hcal_sume_layer_23, "b_hcal_sume_layer_23/F");
    outtree->Branch("hcal_sume_layer_24", &b_hcal_sume_layer_24, "b_hcal_sume_layer_24/F");
    outtree->Branch("hcal_sume_layer_25", &b_hcal_sume_layer_25, "b_hcal_sume_layer_25/F");
    outtree->Branch("hcal_sume_layer_26", &b_hcal_sume_layer_26, "b_hcal_sume_layer_26/F");
    outtree->Branch("hcal_sume_layer_27", &b_hcal_sume_layer_27, "b_hcal_sume_layer_27/F");
    outtree->Branch("hcal_sume_layer_28", &b_hcal_sume_layer_28, "b_hcal_sume_layer_28/F");
    outtree->Branch("hcal_sume_layer_29", &b_hcal_sume_layer_29, "b_hcal_sume_layer_29/F");
    outtree->Branch("hcal_sume_layer_30", &b_hcal_sume_layer_30, "b_hcal_sume_layer_30/F");
    outtree->Branch("hcal_sume_layer_31", &b_hcal_sume_layer_31, "b_hcal_sume_layer_31/F");
    outtree->Branch("hcal_sume_layer_32", &b_hcal_sume_layer_32, "b_hcal_sume_layer_32/F");
    outtree->Branch("hcal_sume_layer_33", &b_hcal_sume_layer_33, "b_hcal_sume_layer_33/F");
    outtree->Branch("hcal_sume_layer_34", &b_hcal_sume_layer_34, "b_hcal_sume_layer_34/F");
    outtree->Branch("hcal_sume_layer_35", &b_hcal_sume_layer_35, "b_hcal_sume_layer_35/F");
    outtree->Branch("hcal_sume_layer_36", &b_hcal_sume_layer_36, "b_hcal_sume_layer_36/F");
    outtree->Branch("hcal_sume_layer_37", &b_hcal_sume_layer_37, "b_hcal_sume_layer_37/F");
    outtree->Branch("hcal_sume_layer_38", &b_hcal_sume_layer_38, "b_hcal_sume_layer_38/F");
    outtree->Branch("hcal_sume_layer_39", &b_hcal_sume_layer_39, "b_hcal_sume_layer_39/F");
    outtree->Branch("hcal_sume_layer_40", &b_hcal_sume_layer_40, "b_hcal_sume_layer_40/F");
    outtree->Branch("hcal_sume_layer_n_0", &b_hcal_sume_layer_n_0, "b_hcal_sume_layer_n_0/F");
    outtree->Branch("hcal_sume_layer_n_1", &b_hcal_sume_layer_n_1, "b_hcal_sume_layer_n_1/F");
    outtree->Branch("hcal_sume_layer_n_2", &b_hcal_sume_layer_n_2, "b_hcal_sume_layer_n_2/F");
    outtree->Branch("hcal_sume_layer_n_3", &b_hcal_sume_layer_n_3, "b_hcal_sume_layer_n_3/F");
    outtree->Branch("hcal_sume_layer_n_4", &b_hcal_sume_layer_n_4, "b_hcal_sume_layer_n_4/F");
    outtree->Branch("hcal_sume_layer_n_5", &b_hcal_sume_layer_n_5, "b_hcal_sume_layer_n_5/F");
    outtree->Branch("hcal_sume_layer_n_6", &b_hcal_sume_layer_n_6, "b_hcal_sume_layer_n_6/F");
    outtree->Branch("hcal_sume_layer_n_7", &b_hcal_sume_layer_n_7, "b_hcal_sume_layer_n_7/F");
    outtree->Branch("hcal_sume_layer_n_8", &b_hcal_sume_layer_n_8, "b_hcal_sume_layer_n_8/F");
    outtree->Branch("hcal_sume_layer_n_9", &b_hcal_sume_layer_n_9, "b_hcal_sume_layer_n_9/F");
    outtree->Branch("hcal_sume_layer_n_10", &b_hcal_sume_layer_n_10, "b_hcal_sume_layer_n_10/F");
    outtree->Branch("hcal_sume_layer_n_11", &b_hcal_sume_layer_n_11, "b_hcal_sume_layer_n_11/F");
    outtree->Branch("hcal_sume_layer_n_12", &b_hcal_sume_layer_n_12, "b_hcal_sume_layer_n_12/F");
    outtree->Branch("hcal_sume_layer_n_13", &b_hcal_sume_layer_n_13, "b_hcal_sume_layer_n_13/F");
    outtree->Branch("hcal_sume_layer_n_14", &b_hcal_sume_layer_n_14, "b_hcal_sume_layer_n_14/F");
    outtree->Branch("hcal_sume_layer_n_15", &b_hcal_sume_layer_n_15, "b_hcal_sume_layer_n_15/F");
    outtree->Branch("hcal_sume_layer_n_16", &b_hcal_sume_layer_n_16, "b_hcal_sume_layer_n_16/F");
    outtree->Branch("hcal_sume_layer_n_17", &b_hcal_sume_layer_n_17, "b_hcal_sume_layer_n_17/F");
    outtree->Branch("hcal_sume_layer_n_18", &b_hcal_sume_layer_n_18, "b_hcal_sume_layer_n_18/F");
    outtree->Branch("hcal_sume_layer_n_19", &b_hcal_sume_layer_n_19, "b_hcal_sume_layer_n_19/F");
    outtree->Branch("hcal_sume_layer_n_20", &b_hcal_sume_layer_n_20, "b_hcal_sume_layer_n_20/F");
    outtree->Branch("hcal_sume_layer_n_21", &b_hcal_sume_layer_n_21, "b_hcal_sume_layer_n_21/F");
    outtree->Branch("hcal_sume_layer_n_22", &b_hcal_sume_layer_n_22, "b_hcal_sume_layer_n_22/F");
    outtree->Branch("hcal_sume_layer_n_23", &b_hcal_sume_layer_n_23, "b_hcal_sume_layer_n_23/F");
    outtree->Branch("hcal_sume_layer_n_24", &b_hcal_sume_layer_n_24, "b_hcal_sume_layer_n_24/F");
    outtree->Branch("hcal_sume_layer_n_25", &b_hcal_sume_layer_n_25, "b_hcal_sume_layer_n_25/F");
    outtree->Branch("hcal_sume_layer_n_26", &b_hcal_sume_layer_n_26, "b_hcal_sume_layer_n_26/F");
    outtree->Branch("hcal_sume_layer_n_27", &b_hcal_sume_layer_n_27, "b_hcal_sume_layer_n_27/F");
    outtree->Branch("hcal_sume_layer_n_28", &b_hcal_sume_layer_n_28, "b_hcal_sume_layer_n_28/F");
    outtree->Branch("hcal_sume_layer_n_29", &b_hcal_sume_layer_n_29, "b_hcal_sume_layer_n_29/F");
    outtree->Branch("hcal_sume_layer_n_30", &b_hcal_sume_layer_n_30, "b_hcal_sume_layer_n_30/F");
    outtree->Branch("hcal_sume_layer_n_31", &b_hcal_sume_layer_n_31, "b_hcal_sume_layer_n_31/F");
    outtree->Branch("hcal_sume_layer_n_32", &b_hcal_sume_layer_n_32, "b_hcal_sume_layer_n_32/F");
    outtree->Branch("hcal_sume_layer_n_33", &b_hcal_sume_layer_n_33, "b_hcal_sume_layer_n_33/F");
    outtree->Branch("hcal_sume_layer_n_34", &b_hcal_sume_layer_n_34, "b_hcal_sume_layer_n_34/F");
    outtree->Branch("hcal_sume_layer_n_35", &b_hcal_sume_layer_n_35, "b_hcal_sume_layer_n_35/F");
    outtree->Branch("hcal_sume_layer_n_36", &b_hcal_sume_layer_n_36, "b_hcal_sume_layer_n_36/F");
    outtree->Branch("hcal_sume_layer_n_37", &b_hcal_sume_layer_n_37, "b_hcal_sume_layer_n_37/F");
    outtree->Branch("hcal_sume_layer_n_38", &b_hcal_sume_layer_n_38, "b_hcal_sume_layer_n_38/F");
    outtree->Branch("hcal_sume_layer_n_39", &b_hcal_sume_layer_n_39, "b_hcal_sume_layer_n_39/F");
    outtree->Branch("hcal_sume_layer_n_40", &b_hcal_sume_layer_n_40, "b_hcal_sume_layer_n_40/F");
    outtree->Branch("hcal_weighte_layer_0", &b_hcal_weighte_layer_0, "b_hcal_weighte_layer_0/F");
    outtree->Branch("hcal_weighte_layer_1", &b_hcal_weighte_layer_1, "b_hcal_weighte_layer_1/F");
    outtree->Branch("hcal_weighte_layer_2", &b_hcal_weighte_layer_2, "b_hcal_weighte_layer_2/F");
    outtree->Branch("hcal_weighte_layer_3", &b_hcal_weighte_layer_3, "b_hcal_weighte_layer_3/F");
    outtree->Branch("hcal_weighte_layer_4", &b_hcal_weighte_layer_4, "b_hcal_weighte_layer_4/F");
    outtree->Branch("hcal_weighte_layer_5", &b_hcal_weighte_layer_5, "b_hcal_weighte_layer_5/F");
    outtree->Branch("hcal_weighte_layer_6", &b_hcal_weighte_layer_6, "b_hcal_weighte_layer_6/F");
    outtree->Branch("hcal_weighte_layer_7", &b_hcal_weighte_layer_7, "b_hcal_weighte_layer_7/F");
    outtree->Branch("hcal_weighte_layer_8", &b_hcal_weighte_layer_8, "b_hcal_weighte_layer_8/F");
    outtree->Branch("hcal_weighte_layer_9", &b_hcal_weighte_layer_9, "b_hcal_weighte_layer_9/F");
    outtree->Branch("hcal_weighte_layer_10", &b_hcal_weighte_layer_10, "b_hcal_weighte_layer_10/F");
    outtree->Branch("hcal_weighte_layer_11", &b_hcal_weighte_layer_11, "b_hcal_weighte_layer_11/F");
    outtree->Branch("hcal_weighte_layer_12", &b_hcal_weighte_layer_12, "b_hcal_weighte_layer_12/F");
    outtree->Branch("hcal_weighte_layer_13", &b_hcal_weighte_layer_13, "b_hcal_weighte_layer_13/F");
    outtree->Branch("hcal_weighte_layer_14", &b_hcal_weighte_layer_14, "b_hcal_weighte_layer_14/F");
    outtree->Branch("hcal_weighte_layer_15", &b_hcal_weighte_layer_15, "b_hcal_weighte_layer_15/F");
    outtree->Branch("hcal_weighte_layer_16", &b_hcal_weighte_layer_16, "b_hcal_weighte_layer_16/F");
    outtree->Branch("hcal_weighte_layer_17", &b_hcal_weighte_layer_17, "b_hcal_weighte_layer_17/F");
    outtree->Branch("hcal_weighte_layer_18", &b_hcal_weighte_layer_18, "b_hcal_weighte_layer_18/F");
    outtree->Branch("hcal_weighte_layer_19", &b_hcal_weighte_layer_19, "b_hcal_weighte_layer_19/F");
    outtree->Branch("hcal_weighte_layer_20", &b_hcal_weighte_layer_20, "b_hcal_weighte_layer_20/F");
    outtree->Branch("hcal_weighte_layer_21", &b_hcal_weighte_layer_21, "b_hcal_weighte_layer_21/F");
    outtree->Branch("hcal_weighte_layer_22", &b_hcal_weighte_layer_22, "b_hcal_weighte_layer_22/F");
    outtree->Branch("hcal_weighte_layer_23", &b_hcal_weighte_layer_23, "b_hcal_weighte_layer_23/F");
    outtree->Branch("hcal_weighte_layer_24", &b_hcal_weighte_layer_24, "b_hcal_weighte_layer_24/F");
    outtree->Branch("hcal_weighte_layer_25", &b_hcal_weighte_layer_25, "b_hcal_weighte_layer_25/F");
    outtree->Branch("hcal_weighte_layer_26", &b_hcal_weighte_layer_26, "b_hcal_weighte_layer_26/F");
    outtree->Branch("hcal_weighte_layer_27", &b_hcal_weighte_layer_27, "b_hcal_weighte_layer_27/F");
    outtree->Branch("hcal_weighte_layer_28", &b_hcal_weighte_layer_28, "b_hcal_weighte_layer_28/F");
    outtree->Branch("hcal_weighte_layer_29", &b_hcal_weighte_layer_29, "b_hcal_weighte_layer_29/F");
    outtree->Branch("hcal_weighte_layer_30", &b_hcal_weighte_layer_30, "b_hcal_weighte_layer_30/F");
    outtree->Branch("hcal_weighte_layer_31", &b_hcal_weighte_layer_31, "b_hcal_weighte_layer_31/F");
    outtree->Branch("hcal_weighte_layer_32", &b_hcal_weighte_layer_32, "b_hcal_weighte_layer_32/F");
    outtree->Branch("hcal_weighte_layer_33", &b_hcal_weighte_layer_33, "b_hcal_weighte_layer_33/F");
    outtree->Branch("hcal_weighte_layer_34", &b_hcal_weighte_layer_34, "b_hcal_weighte_layer_34/F");
    outtree->Branch("hcal_weighte_layer_35", &b_hcal_weighte_layer_35, "b_hcal_weighte_layer_35/F");
    outtree->Branch("hcal_weighte_layer_36", &b_hcal_weighte_layer_36, "b_hcal_weighte_layer_36/F");
    outtree->Branch("hcal_weighte_layer_37", &b_hcal_weighte_layer_37, "b_hcal_weighte_layer_37/F");
    outtree->Branch("hcal_weighte_layer_38", &b_hcal_weighte_layer_38, "b_hcal_weighte_layer_38/F");
    outtree->Branch("hcal_weighte_layer_39", &b_hcal_weighte_layer_39, "b_hcal_weighte_layer_39/F");
    outtree->Branch("hcal_weighte_layer_40", &b_hcal_weighte_layer_40, "b_hcal_weighte_layer_40/F");
    outtree->Branch("hcal_weighte_layer_n_0", &b_hcal_weighte_layer_n_0, "b_hcal_weighte_layer_n_0/F");
    outtree->Branch("hcal_weighte_layer_n_1", &b_hcal_weighte_layer_n_1, "b_hcal_weighte_layer_n_1/F");
    outtree->Branch("hcal_weighte_layer_n_2", &b_hcal_weighte_layer_n_2, "b_hcal_weighte_layer_n_2/F");
    outtree->Branch("hcal_weighte_layer_n_3", &b_hcal_weighte_layer_n_3, "b_hcal_weighte_layer_n_3/F");
    outtree->Branch("hcal_weighte_layer_n_4", &b_hcal_weighte_layer_n_4, "b_hcal_weighte_layer_n_4/F");
    outtree->Branch("hcal_weighte_layer_n_5", &b_hcal_weighte_layer_n_5, "b_hcal_weighte_layer_n_5/F");
    outtree->Branch("hcal_weighte_layer_n_6", &b_hcal_weighte_layer_n_6, "b_hcal_weighte_layer_n_6/F");
    outtree->Branch("hcal_weighte_layer_n_7", &b_hcal_weighte_layer_n_7, "b_hcal_weighte_layer_n_7/F");
    outtree->Branch("hcal_weighte_layer_n_8", &b_hcal_weighte_layer_n_8, "b_hcal_weighte_layer_n_8/F");
    outtree->Branch("hcal_weighte_layer_n_9", &b_hcal_weighte_layer_n_9, "b_hcal_weighte_layer_n_9/F");
    outtree->Branch("hcal_weighte_layer_n_10", &b_hcal_weighte_layer_n_10, "b_hcal_weighte_layer_n_10/F");
    outtree->Branch("hcal_weighte_layer_n_11", &b_hcal_weighte_layer_n_11, "b_hcal_weighte_layer_n_11/F");
    outtree->Branch("hcal_weighte_layer_n_12", &b_hcal_weighte_layer_n_12, "b_hcal_weighte_layer_n_12/F");
    outtree->Branch("hcal_weighte_layer_n_13", &b_hcal_weighte_layer_n_13, "b_hcal_weighte_layer_n_13/F");
    outtree->Branch("hcal_weighte_layer_n_14", &b_hcal_weighte_layer_n_14, "b_hcal_weighte_layer_n_14/F");
    outtree->Branch("hcal_weighte_layer_n_15", &b_hcal_weighte_layer_n_15, "b_hcal_weighte_layer_n_15/F");
    outtree->Branch("hcal_weighte_layer_n_16", &b_hcal_weighte_layer_n_16, "b_hcal_weighte_layer_n_16/F");
    outtree->Branch("hcal_weighte_layer_n_17", &b_hcal_weighte_layer_n_17, "b_hcal_weighte_layer_n_17/F");
    outtree->Branch("hcal_weighte_layer_n_18", &b_hcal_weighte_layer_n_18, "b_hcal_weighte_layer_n_18/F");
    outtree->Branch("hcal_weighte_layer_n_19", &b_hcal_weighte_layer_n_19, "b_hcal_weighte_layer_n_19/F");
    outtree->Branch("hcal_weighte_layer_n_20", &b_hcal_weighte_layer_n_20, "b_hcal_weighte_layer_n_20/F");
    outtree->Branch("hcal_weighte_layer_n_21", &b_hcal_weighte_layer_n_21, "b_hcal_weighte_layer_n_21/F");
    outtree->Branch("hcal_weighte_layer_n_22", &b_hcal_weighte_layer_n_22, "b_hcal_weighte_layer_n_22/F");
    outtree->Branch("hcal_weighte_layer_n_23", &b_hcal_weighte_layer_n_23, "b_hcal_weighte_layer_n_23/F");
    outtree->Branch("hcal_weighte_layer_n_24", &b_hcal_weighte_layer_n_24, "b_hcal_weighte_layer_n_24/F");
    outtree->Branch("hcal_weighte_layer_n_25", &b_hcal_weighte_layer_n_25, "b_hcal_weighte_layer_n_25/F");
    outtree->Branch("hcal_weighte_layer_n_26", &b_hcal_weighte_layer_n_26, "b_hcal_weighte_layer_n_26/F");
    outtree->Branch("hcal_weighte_layer_n_27", &b_hcal_weighte_layer_n_27, "b_hcal_weighte_layer_n_27/F");
    outtree->Branch("hcal_weighte_layer_n_28", &b_hcal_weighte_layer_n_28, "b_hcal_weighte_layer_n_28/F");
    outtree->Branch("hcal_weighte_layer_n_29", &b_hcal_weighte_layer_n_29, "b_hcal_weighte_layer_n_29/F");
    outtree->Branch("hcal_weighte_layer_n_30", &b_hcal_weighte_layer_n_30, "b_hcal_weighte_layer_n_30/F");
    outtree->Branch("hcal_weighte_layer_n_31", &b_hcal_weighte_layer_n_31, "b_hcal_weighte_layer_n_31/F");
    outtree->Branch("hcal_weighte_layer_n_32", &b_hcal_weighte_layer_n_32, "b_hcal_weighte_layer_n_32/F");
    outtree->Branch("hcal_weighte_layer_n_33", &b_hcal_weighte_layer_n_33, "b_hcal_weighte_layer_n_33/F");
    outtree->Branch("hcal_weighte_layer_n_34", &b_hcal_weighte_layer_n_34, "b_hcal_weighte_layer_n_34/F");
    outtree->Branch("hcal_weighte_layer_n_35", &b_hcal_weighte_layer_n_35, "b_hcal_weighte_layer_n_35/F");
    outtree->Branch("hcal_weighte_layer_n_36", &b_hcal_weighte_layer_n_36, "b_hcal_weighte_layer_n_36/F");
    outtree->Branch("hcal_weighte_layer_n_37", &b_hcal_weighte_layer_n_37, "b_hcal_weighte_layer_n_37/F");
    outtree->Branch("hcal_weighte_layer_n_38", &b_hcal_weighte_layer_n_38, "b_hcal_weighte_layer_n_38/F");
    outtree->Branch("hcal_weighte_layer_n_39", &b_hcal_weighte_layer_n_39, "b_hcal_weighte_layer_n_39/F");
    outtree->Branch("hcal_weighte_layer_n_40", &b_hcal_weighte_layer_n_40, "b_hcal_weighte_layer_n_40/F");

    outtree->Branch("total_interaction",&b_total_interaction,"b_total_interaction/F");
    outtree->Branch("total_nhit",&b_total_nhit,"b_total_nhit/F");
    outtree->Branch("total_sume",&b_total_sume,"b_total_sume/F");
    outtree->Branch("total_weighte",&b_total_weighte,"b_total_weighte/F");
    outtree->Branch("total_mol",&b_total_mol,"b_total_mol/F");
    outtree->Branch("total_MIP_Likeness",&b_total_MIP_Likeness,"b_total_MIP_Likeness/F");
    outtree->Branch("total_hits_max_distance",&b_total_hits_max_distance,"b_total_hits_max_distance/F");
    outtree->Branch("total_radius90_layer_0",&b_total_radius90_layer_0,"b_total_radius90_layer_0/F");
    outtree->Branch("total_radius90_layer_1",&b_total_radius90_layer_1,"b_total_radius90_layer_1/F");
    outtree->Branch("total_radius90_layer_2",&b_total_radius90_layer_2,"b_total_radius90_layer_2/F");
    outtree->Branch("total_radius90_layer_3",&b_total_radius90_layer_3,"b_total_radius90_layer_3/F");
    outtree->Branch("total_radius90_layer_4",&b_total_radius90_layer_4,"b_total_radius90_layer_4/F");
    outtree->Branch("total_radius90_layer_5",&b_total_radius90_layer_5,"b_total_radius90_layer_5/F");
    outtree->Branch("total_radius90_layer_6",&b_total_radius90_layer_6,"b_total_radius90_layer_6/F");
    outtree->Branch("total_radius90_layer_7",&b_total_radius90_layer_7,"b_total_radius90_layer_7/F");
    outtree->Branch("total_radius90_layer_8",&b_total_radius90_layer_8,"b_total_radius90_layer_8/F");
    outtree->Branch("total_radius90_layer_9",&b_total_radius90_layer_9,"b_total_radius90_layer_9/F");
    outtree->Branch("total_radius90_layer_10",&b_total_radius90_layer_10,"b_total_radius90_layer_10/F");
    outtree->Branch("total_radius90_layer_11",&b_total_radius90_layer_11,"b_total_radius90_layer_11/F");
    outtree->Branch("total_radius90_layer_12",&b_total_radius90_layer_12,"b_total_radius90_layer_12/F");
    outtree->Branch("total_radius90_layer_13",&b_total_radius90_layer_13,"b_total_radius90_layer_13/F");
    outtree->Branch("total_radius90_layer_14",&b_total_radius90_layer_14,"b_total_radius90_layer_14/F");
    outtree->Branch("total_bar_x",&b_total_bar_x,"b_total_bar_x/F");
    outtree->Branch("total_bar_y",&b_total_bar_y,"b_total_bar_y/F");
    outtree->Branch("total_bar_z",&b_total_bar_z,"b_total_bar_z/F");
    outtree->Branch("total_bar_r",&b_total_bar_r,"b_total_bar_r/F");
    outtree->Branch("total_bar_x_layer_0",&b_total_bar_x_layer_0,"b_total_bar_x_layer_0/F");
    outtree->Branch("total_bar_x_layer_1",&b_total_bar_x_layer_1,"b_total_bar_x_layer_1/F");
    outtree->Branch("total_bar_x_layer_2",&b_total_bar_x_layer_2,"b_total_bar_x_layer_2/F");
    outtree->Branch("total_bar_x_layer_3",&b_total_bar_x_layer_3,"b_total_bar_x_layer_3/F");
    outtree->Branch("total_bar_x_layer_4",&b_total_bar_x_layer_4,"b_total_bar_x_layer_4/F");
    outtree->Branch("total_bar_x_layer_5",&b_total_bar_x_layer_5,"b_total_bar_x_layer_5/F");
    outtree->Branch("total_bar_x_layer_6",&b_total_bar_x_layer_6,"b_total_bar_x_layer_6/F");
    outtree->Branch("total_bar_x_layer_7",&b_total_bar_x_layer_7,"b_total_bar_x_layer_7/F");
    outtree->Branch("total_bar_x_layer_8",&b_total_bar_x_layer_8,"b_total_bar_x_layer_8/F");
    outtree->Branch("total_bar_x_layer_9",&b_total_bar_x_layer_9,"b_total_bar_x_layer_9/F");
    outtree->Branch("total_bar_x_layer_10",&b_total_bar_x_layer_10,"b_total_bar_x_layer_10/F");
    outtree->Branch("total_bar_x_layer_11",&b_total_bar_x_layer_11,"b_total_bar_x_layer_11/F");
    outtree->Branch("total_bar_x_layer_12",&b_total_bar_x_layer_12,"b_total_bar_x_layer_12/F");
    outtree->Branch("total_bar_x_layer_13",&b_total_bar_x_layer_13,"b_total_bar_x_layer_13/F");
    outtree->Branch("total_bar_x_layer_14",&b_total_bar_x_layer_14,"b_total_bar_x_layer_14/F");
    outtree->Branch("total_bar_y_layer_0",&b_total_bar_y_layer_0,"b_total_bar_y_layer_0/F");
    outtree->Branch("total_bar_y_layer_1",&b_total_bar_y_layer_1,"b_total_bar_y_layer_1/F");
    outtree->Branch("total_bar_y_layer_2",&b_total_bar_y_layer_2,"b_total_bar_y_layer_2/F");
    outtree->Branch("total_bar_y_layer_3",&b_total_bar_y_layer_3,"b_total_bar_y_layer_3/F");
    outtree->Branch("total_bar_y_layer_4",&b_total_bar_y_layer_4,"b_total_bar_y_layer_4/F");
    outtree->Branch("total_bar_y_layer_5",&b_total_bar_y_layer_5,"b_total_bar_y_layer_5/F");
    outtree->Branch("total_bar_y_layer_6",&b_total_bar_y_layer_6,"b_total_bar_y_layer_6/F");
    outtree->Branch("total_bar_y_layer_7",&b_total_bar_y_layer_7,"b_total_bar_y_layer_7/F");
    outtree->Branch("total_bar_y_layer_8",&b_total_bar_y_layer_8,"b_total_bar_y_layer_8/F");
    outtree->Branch("total_bar_y_layer_9",&b_total_bar_y_layer_9,"b_total_bar_y_layer_9/F");
    outtree->Branch("total_bar_y_layer_10",&b_total_bar_y_layer_10,"b_total_bar_y_layer_10/F");
    outtree->Branch("total_bar_y_layer_11",&b_total_bar_y_layer_11,"b_total_bar_y_layer_11/F");
    outtree->Branch("total_bar_y_layer_12",&b_total_bar_y_layer_12,"b_total_bar_y_layer_12/F");
    outtree->Branch("total_bar_y_layer_13",&b_total_bar_y_layer_13,"b_total_bar_y_layer_13/F");
    outtree->Branch("total_bar_y_layer_14",&b_total_bar_y_layer_14,"b_total_bar_y_layer_14/F");
    outtree->Branch("total_bar_r_layer_0",&b_total_bar_r_layer_0,"b_total_bar_r_layer_0/F");
    outtree->Branch("total_bar_r_layer_1",&b_total_bar_r_layer_1,"b_total_bar_r_layer_1/F");
    outtree->Branch("total_bar_r_layer_2",&b_total_bar_r_layer_2,"b_total_bar_r_layer_2/F");
    outtree->Branch("total_bar_r_layer_3",&b_total_bar_r_layer_3,"b_total_bar_r_layer_3/F");
    outtree->Branch("total_bar_r_layer_4",&b_total_bar_r_layer_4,"b_total_bar_r_layer_4/F");
    outtree->Branch("total_bar_r_layer_5",&b_total_bar_r_layer_5,"b_total_bar_r_layer_5/F");
    outtree->Branch("total_bar_r_layer_6",&b_total_bar_r_layer_6,"b_total_bar_r_layer_6/F");
    outtree->Branch("total_bar_r_layer_7",&b_total_bar_r_layer_7,"b_total_bar_r_layer_7/F");
    outtree->Branch("total_bar_r_layer_8",&b_total_bar_r_layer_8,"b_total_bar_r_layer_8/F");
    outtree->Branch("total_bar_r_layer_9",&b_total_bar_r_layer_9,"b_total_bar_r_layer_9/F");
    outtree->Branch("total_bar_r_layer_10",&b_total_bar_r_layer_10,"b_total_bar_r_layer_10/F");
    outtree->Branch("total_bar_r_layer_11",&b_total_bar_r_layer_11,"b_total_bar_r_layer_11/F");
    outtree->Branch("total_bar_r_layer_12",&b_total_bar_r_layer_12,"b_total_bar_r_layer_12/F");
    outtree->Branch("total_bar_r_layer_13",&b_total_bar_r_layer_13,"b_total_bar_r_layer_13/F");
    outtree->Branch("total_bar_r_layer_14",&b_total_bar_r_layer_14,"b_total_bar_r_layer_14/F");
    outtree->Branch("total_shower_nhit_max_layer",&b_total_shower_nhit_max_layer,"b_total_shower_nhit_max_layer/F");
    outtree->Branch("total_shower_nhit_start_layer",&b_total_shower_nhit_start_layer,"b_total_shower_nhit_start_layer/F");
    outtree->Branch("total_shower_nhit_end_layer",&b_total_shower_nhit_end_layer,"b_total_shower_nhit_end_layer/F");
    outtree->Branch("total_shower_nhit_start_10_layer",&b_total_shower_nhit_start_10_layer,"b_total_shower_nhit_start_10_layer/F");
    outtree->Branch("total_shower_nhit_end_10_layer",&b_total_shower_nhit_end_10_layer,"b_total_shower_nhit_end_10_layer/F");
    outtree->Branch("total_shower_nhit_average",&b_total_shower_nhit_average,"b_total_shower_nhit_average/F");
    outtree->Branch("total_shower_nhit_max",&b_total_shower_nhit_max,"b_total_shower_nhit_max/F");
    outtree->Branch("total_shower_sume_max_layer",&b_total_shower_sume_max_layer,"b_total_shower_sume_max_layer/F");
    outtree->Branch("total_shower_sume_start_layer",&b_total_shower_sume_start_layer,"b_total_shower_sume_start_layer/F");
    outtree->Branch("total_shower_sume_end_layer",&b_total_shower_sume_end_layer,"b_total_shower_sume_end_layer/F");
    outtree->Branch("total_shower_sume_start_10_layer",&b_total_shower_sume_start_10_layer,"b_total_shower_sume_start_10_layer/F");
    outtree->Branch("total_shower_sume_end_10_layer",&b_total_shower_sume_end_10_layer,"b_total_shower_sume_end_10_layer/F");
    outtree->Branch("total_shower_sume_average",&b_total_shower_sume_average,"b_total_shower_sume_average/F");
    outtree->Branch("total_shower_sume_max",&b_total_shower_sume_max,"b_total_shower_sume_max/F");
    outtree->Branch("total_shower_weighte_max_layer",&b_total_shower_weighte_max_layer,"b_total_shower_weighte_max_layer/F");
    outtree->Branch("total_shower_weighte_start_layer",&b_total_shower_weighte_start_layer,"b_total_shower_weighte_start_layer/F");
    outtree->Branch("total_shower_weighte_end_layer",&b_total_shower_weighte_end_layer,"b_total_shower_weighte_end_layer/F");
    outtree->Branch("total_shower_weighte_start_10_layer",&b_total_shower_weighte_start_10_layer,"b_total_shower_weighte_start_10_layer/F");
    outtree->Branch("total_shower_weighte_end_10_layer",&b_total_shower_weighte_end_10_layer,"b_total_shower_weighte_end_10_layer/F");
    outtree->Branch("total_shower_weighte_average",&b_total_shower_weighte_average,"b_total_shower_weighte_average/F");
    outtree->Branch("total_shower_weighte_max",&b_total_shower_weighte_max,"b_total_shower_weighte_max/F");
    outtree->Branch("total_nhit_layer_0", &b_total_nhit_layer_0, "b_total_nhit_layer_0/F");
    outtree->Branch("total_nhit_layer_1", &b_total_nhit_layer_1, "b_total_nhit_layer_1/F");
    outtree->Branch("total_nhit_layer_2", &b_total_nhit_layer_2, "b_total_nhit_layer_2/F");
    outtree->Branch("total_nhit_layer_3", &b_total_nhit_layer_3, "b_total_nhit_layer_3/F");
    outtree->Branch("total_nhit_layer_4", &b_total_nhit_layer_4, "b_total_nhit_layer_4/F");
    outtree->Branch("total_nhit_layer_5", &b_total_nhit_layer_5, "b_total_nhit_layer_5/F");
    outtree->Branch("total_nhit_layer_6", &b_total_nhit_layer_6, "b_total_nhit_layer_6/F");
    outtree->Branch("total_nhit_layer_7", &b_total_nhit_layer_7, "b_total_nhit_layer_7/F");
    outtree->Branch("total_nhit_layer_8", &b_total_nhit_layer_8, "b_total_nhit_layer_8/F");
    outtree->Branch("total_nhit_layer_9", &b_total_nhit_layer_9, "b_total_nhit_layer_9/F");
    outtree->Branch("total_nhit_layer_10", &b_total_nhit_layer_10, "b_total_nhit_layer_10/F");
    outtree->Branch("total_nhit_layer_11", &b_total_nhit_layer_11, "b_total_nhit_layer_11/F");
    outtree->Branch("total_nhit_layer_12", &b_total_nhit_layer_12, "b_total_nhit_layer_12/F");
    outtree->Branch("total_nhit_layer_13", &b_total_nhit_layer_13, "b_total_nhit_layer_13/F");
    outtree->Branch("total_nhit_layer_14", &b_total_nhit_layer_14, "b_total_nhit_layer_14/F");
    outtree->Branch("total_nhit_layer_15", &b_total_nhit_layer_15, "b_total_nhit_layer_15/F");
    outtree->Branch("total_nhit_layer_16", &b_total_nhit_layer_16, "b_total_nhit_layer_16/F");
    outtree->Branch("total_nhit_layer_17", &b_total_nhit_layer_17, "b_total_nhit_layer_17/F");
    outtree->Branch("total_nhit_layer_18", &b_total_nhit_layer_18, "b_total_nhit_layer_18/F");
    outtree->Branch("total_nhit_layer_19", &b_total_nhit_layer_19, "b_total_nhit_layer_19/F");
    outtree->Branch("total_nhit_layer_20", &b_total_nhit_layer_20, "b_total_nhit_layer_20/F");
    outtree->Branch("total_nhit_layer_21", &b_total_nhit_layer_21, "b_total_nhit_layer_21/F");
    outtree->Branch("total_nhit_layer_22", &b_total_nhit_layer_22, "b_total_nhit_layer_22/F");
    outtree->Branch("total_nhit_layer_23", &b_total_nhit_layer_23, "b_total_nhit_layer_23/F");
    outtree->Branch("total_nhit_layer_24", &b_total_nhit_layer_24, "b_total_nhit_layer_24/F");
    outtree->Branch("total_nhit_layer_25", &b_total_nhit_layer_25, "b_total_nhit_layer_25/F");
    outtree->Branch("total_nhit_layer_26", &b_total_nhit_layer_26, "b_total_nhit_layer_26/F");
    outtree->Branch("total_nhit_layer_27", &b_total_nhit_layer_27, "b_total_nhit_layer_27/F");
    outtree->Branch("total_nhit_layer_28", &b_total_nhit_layer_28, "b_total_nhit_layer_28/F");
    outtree->Branch("total_nhit_layer_29", &b_total_nhit_layer_29, "b_total_nhit_layer_29/F");
    outtree->Branch("total_nhit_layer_30", &b_total_nhit_layer_30, "b_total_nhit_layer_30/F");
    outtree->Branch("total_nhit_layer_31", &b_total_nhit_layer_31, "b_total_nhit_layer_31/F");
    outtree->Branch("total_nhit_layer_32", &b_total_nhit_layer_32, "b_total_nhit_layer_32/F");
    outtree->Branch("total_nhit_layer_33", &b_total_nhit_layer_33, "b_total_nhit_layer_33/F");
    outtree->Branch("total_nhit_layer_34", &b_total_nhit_layer_34, "b_total_nhit_layer_34/F");
    outtree->Branch("total_nhit_layer_35", &b_total_nhit_layer_35, "b_total_nhit_layer_35/F");
    outtree->Branch("total_nhit_layer_36", &b_total_nhit_layer_36, "b_total_nhit_layer_36/F");
    outtree->Branch("total_nhit_layer_37", &b_total_nhit_layer_37, "b_total_nhit_layer_37/F");
    outtree->Branch("total_nhit_layer_38", &b_total_nhit_layer_38, "b_total_nhit_layer_38/F");
    outtree->Branch("total_nhit_layer_39", &b_total_nhit_layer_39, "b_total_nhit_layer_39/F");
    outtree->Branch("total_nhit_layer_40", &b_total_nhit_layer_40, "b_total_nhit_layer_40/F");
    outtree->Branch("total_nhit_layer_41", &b_total_nhit_layer_41, "b_total_nhit_layer_41/F");
    outtree->Branch("total_nhit_layer_42", &b_total_nhit_layer_42, "b_total_nhit_layer_42/F");
    outtree->Branch("total_nhit_layer_43", &b_total_nhit_layer_43, "b_total_nhit_layer_43F");
    outtree->Branch("total_nhit_layer_44", &b_total_nhit_layer_44, "b_total_nhit_layer_44F");
    outtree->Branch("total_nhit_layer_45", &b_total_nhit_layer_45, "b_total_nhit_layer_45F");
    outtree->Branch("total_nhit_layer_46", &b_total_nhit_layer_46, "b_total_nhit_layer_46F");
    outtree->Branch("total_nhit_layer_47", &b_total_nhit_layer_47, "b_total_nhit_layer_47F");
    outtree->Branch("total_nhit_layer_48", &b_total_nhit_layer_48, "b_total_nhit_layer_48F");
    outtree->Branch("total_nhit_layer_49", &b_total_nhit_layer_49, "b_total_nhit_layer_49F");
    outtree->Branch("total_nhit_layer_50", &b_total_nhit_layer_50, "b_total_nhit_layer_50/F");
    outtree->Branch("total_nhit_layer_51", &b_total_nhit_layer_51, "b_total_nhit_layer_51/F");
    outtree->Branch("total_nhit_layer_52", &b_total_nhit_layer_52, "b_total_nhit_layer_52/F");
    outtree->Branch("total_nhit_layer_53", &b_total_nhit_layer_53, "b_total_nhit_layer_53/F");
    outtree->Branch("total_nhit_layer_54", &b_total_nhit_layer_54, "b_total_nhit_layer_54/F");
    outtree->Branch("total_nhit_layer_55", &b_total_nhit_layer_55, "b_total_nhit_layer_55/F");
    outtree->Branch("total_nhit_layer_n_0", &b_total_nhit_layer_n_0, "b_total_nhit_layer_n_0/F");
    outtree->Branch("total_nhit_layer_n_1", &b_total_nhit_layer_n_1, "b_total_nhit_layer_n_1/F");
    outtree->Branch("total_nhit_layer_n_2", &b_total_nhit_layer_n_2, "b_total_nhit_layer_n_2/F");
    outtree->Branch("total_nhit_layer_n_3", &b_total_nhit_layer_n_3, "b_total_nhit_layer_n_3/F");
    outtree->Branch("total_nhit_layer_n_4", &b_total_nhit_layer_n_4, "b_total_nhit_layer_n_4/F");
    outtree->Branch("total_nhit_layer_n_5", &b_total_nhit_layer_n_5, "b_total_nhit_layer_n_5/F");
    outtree->Branch("total_nhit_layer_n_6", &b_total_nhit_layer_n_6, "b_total_nhit_layer_n_6/F");
    outtree->Branch("total_nhit_layer_n_7", &b_total_nhit_layer_n_7, "b_total_nhit_layer_n_7/F");
    outtree->Branch("total_nhit_layer_n_8", &b_total_nhit_layer_n_8, "b_total_nhit_layer_n_8/F");
    outtree->Branch("total_nhit_layer_n_9", &b_total_nhit_layer_n_9, "b_total_nhit_layer_n_9/F");
    outtree->Branch("total_nhit_layer_n_10", &b_total_nhit_layer_n_10, "b_total_nhit_layer_n_10/F");
    outtree->Branch("total_nhit_layer_n_11", &b_total_nhit_layer_n_11, "b_total_nhit_layer_n_11/F");
    outtree->Branch("total_nhit_layer_n_12", &b_total_nhit_layer_n_12, "b_total_nhit_layer_n_12/F");
    outtree->Branch("total_nhit_layer_n_13", &b_total_nhit_layer_n_13, "b_total_nhit_layer_n_13/F");
    outtree->Branch("total_nhit_layer_n_14", &b_total_nhit_layer_n_14, "b_total_nhit_layer_n_14/F");
    outtree->Branch("total_nhit_layer_n_15", &b_total_nhit_layer_n_15, "b_total_nhit_layer_n_15/F");
    outtree->Branch("total_nhit_layer_n_16", &b_total_nhit_layer_n_16, "b_total_nhit_layer_n_16/F");
    outtree->Branch("total_nhit_layer_n_17", &b_total_nhit_layer_n_17, "b_total_nhit_layer_n_17/F");
    outtree->Branch("total_nhit_layer_n_18", &b_total_nhit_layer_n_18, "b_total_nhit_layer_n_18/F");
    outtree->Branch("total_nhit_layer_n_19", &b_total_nhit_layer_n_19, "b_total_nhit_layer_n_19/F");
    outtree->Branch("total_nhit_layer_n_20", &b_total_nhit_layer_n_20, "b_total_nhit_layer_n_20/F");
    outtree->Branch("total_nhit_layer_n_21", &b_total_nhit_layer_n_21, "b_total_nhit_layer_n_21/F");
    outtree->Branch("total_nhit_layer_n_22", &b_total_nhit_layer_n_22, "b_total_nhit_layer_n_22/F");
    outtree->Branch("total_nhit_layer_n_23", &b_total_nhit_layer_n_23, "b_total_nhit_layer_n_23/F");
    outtree->Branch("total_nhit_layer_n_24", &b_total_nhit_layer_n_24, "b_total_nhit_layer_n_24/F");
    outtree->Branch("total_nhit_layer_n_25", &b_total_nhit_layer_n_25, "b_total_nhit_layer_n_25/F");
    outtree->Branch("total_nhit_layer_n_26", &b_total_nhit_layer_n_26, "b_total_nhit_layer_n_26/F");
    outtree->Branch("total_nhit_layer_n_27", &b_total_nhit_layer_n_27, "b_total_nhit_layer_n_27/F");
    outtree->Branch("total_nhit_layer_n_28", &b_total_nhit_layer_n_28, "b_total_nhit_layer_n_28/F");
    outtree->Branch("total_nhit_layer_n_29", &b_total_nhit_layer_n_29, "b_total_nhit_layer_n_29/F");
    outtree->Branch("total_nhit_layer_n_30", &b_total_nhit_layer_n_30, "b_total_nhit_layer_n_30/F");
    outtree->Branch("total_nhit_layer_n_31", &b_total_nhit_layer_n_31, "b_total_nhit_layer_n_31/F");
    outtree->Branch("total_nhit_layer_n_32", &b_total_nhit_layer_n_32, "b_total_nhit_layer_n_32/F");
    outtree->Branch("total_nhit_layer_n_33", &b_total_nhit_layer_n_33, "b_total_nhit_layer_n_33/F");
    outtree->Branch("total_nhit_layer_n_34", &b_total_nhit_layer_n_34, "b_total_nhit_layer_n_34/F");
    outtree->Branch("total_nhit_layer_n_35", &b_total_nhit_layer_n_35, "b_total_nhit_layer_n_35/F");
    outtree->Branch("total_nhit_layer_n_36", &b_total_nhit_layer_n_36, "b_total_nhit_layer_n_36/F");
    outtree->Branch("total_nhit_layer_n_37", &b_total_nhit_layer_n_37, "b_total_nhit_layer_n_37/F");
    outtree->Branch("total_nhit_layer_n_38", &b_total_nhit_layer_n_38, "b_total_nhit_layer_n_38/F");
    outtree->Branch("total_nhit_layer_n_39", &b_total_nhit_layer_n_39, "b_total_nhit_layer_n_39/F");
    outtree->Branch("total_nhit_layer_n_40", &b_total_nhit_layer_n_40, "b_total_nhit_layer_n_40/F");
    outtree->Branch("total_nhit_layer_n_41", &b_total_nhit_layer_n_41, "b_total_nhit_layer_n_41/F");
    outtree->Branch("total_nhit_layer_n_42", &b_total_nhit_layer_n_42, "b_total_nhit_layer_n_42/F");
    outtree->Branch("total_nhit_layer_n_43", &b_total_nhit_layer_n_43, "b_total_nhit_layer_n_43F");
    outtree->Branch("total_nhit_layer_n_44", &b_total_nhit_layer_n_44, "b_total_nhit_layer_n_44F");
    outtree->Branch("total_nhit_layer_n_45", &b_total_nhit_layer_n_45, "b_total_nhit_layer_n_45F");
    outtree->Branch("total_nhit_layer_n_46", &b_total_nhit_layer_n_46, "b_total_nhit_layer_n_46F");
    outtree->Branch("total_nhit_layer_n_47", &b_total_nhit_layer_n_47, "b_total_nhit_layer_n_47F");
    outtree->Branch("total_nhit_layer_n_48", &b_total_nhit_layer_n_48, "b_total_nhit_layer_n_48F");
    outtree->Branch("total_nhit_layer_n_49", &b_total_nhit_layer_n_49, "b_total_nhit_layer_n_49F");
    outtree->Branch("total_nhit_layer_n_50", &b_total_nhit_layer_n_50, "b_total_nhit_layer_n_50/F");
    outtree->Branch("total_nhit_layer_n_51", &b_total_nhit_layer_n_51, "b_total_nhit_layer_n_51/F");
    outtree->Branch("total_nhit_layer_n_52", &b_total_nhit_layer_n_52, "b_total_nhit_layer_n_52/F");
    outtree->Branch("total_nhit_layer_n_53", &b_total_nhit_layer_n_53, "b_total_nhit_layer_n_53/F");
    outtree->Branch("total_nhit_layer_n_54", &b_total_nhit_layer_n_54, "b_total_nhit_layer_n_54/F");
    outtree->Branch("total_nhit_layer_n_55", &b_total_nhit_layer_n_55, "b_total_nhit_layer_n_55/F");
    outtree->Branch("total_sume_layer_0", &b_total_sume_layer_0, "b_total_sume_layer_0/F");
    outtree->Branch("total_sume_layer_1", &b_total_sume_layer_1, "b_total_sume_layer_1/F");
    outtree->Branch("total_sume_layer_2", &b_total_sume_layer_2, "b_total_sume_layer_2/F");
    outtree->Branch("total_sume_layer_3", &b_total_sume_layer_3, "b_total_sume_layer_3/F");
    outtree->Branch("total_sume_layer_4", &b_total_sume_layer_4, "b_total_sume_layer_4/F");
    outtree->Branch("total_sume_layer_5", &b_total_sume_layer_5, "b_total_sume_layer_5/F");
    outtree->Branch("total_sume_layer_6", &b_total_sume_layer_6, "b_total_sume_layer_6/F");
    outtree->Branch("total_sume_layer_7", &b_total_sume_layer_7, "b_total_sume_layer_7/F");
    outtree->Branch("total_sume_layer_8", &b_total_sume_layer_8, "b_total_sume_layer_8/F");
    outtree->Branch("total_sume_layer_9", &b_total_sume_layer_9, "b_total_sume_layer_9/F");
    outtree->Branch("total_sume_layer_10", &b_total_sume_layer_10, "b_total_sume_layer_10/F");
    outtree->Branch("total_sume_layer_11", &b_total_sume_layer_11, "b_total_sume_layer_11/F");
    outtree->Branch("total_sume_layer_12", &b_total_sume_layer_12, "b_total_sume_layer_12/F");
    outtree->Branch("total_sume_layer_13", &b_total_sume_layer_13, "b_total_sume_layer_13/F");
    outtree->Branch("total_sume_layer_14", &b_total_sume_layer_14, "b_total_sume_layer_14/F");
    outtree->Branch("total_sume_layer_15", &b_total_sume_layer_15, "b_total_sume_layer_15/F");
    outtree->Branch("total_sume_layer_16", &b_total_sume_layer_16, "b_total_sume_layer_16/F");
    outtree->Branch("total_sume_layer_17", &b_total_sume_layer_17, "b_total_sume_layer_17/F");
    outtree->Branch("total_sume_layer_18", &b_total_sume_layer_18, "b_total_sume_layer_18/F");
    outtree->Branch("total_sume_layer_19", &b_total_sume_layer_19, "b_total_sume_layer_19/F");
    outtree->Branch("total_sume_layer_20", &b_total_sume_layer_20, "b_total_sume_layer_20/F");
    outtree->Branch("total_sume_layer_21", &b_total_sume_layer_21, "b_total_sume_layer_21/F");
    outtree->Branch("total_sume_layer_22", &b_total_sume_layer_22, "b_total_sume_layer_22/F");
    outtree->Branch("total_sume_layer_23", &b_total_sume_layer_23, "b_total_sume_layer_23/F");
    outtree->Branch("total_sume_layer_24", &b_total_sume_layer_24, "b_total_sume_layer_24/F");
    outtree->Branch("total_sume_layer_25", &b_total_sume_layer_25, "b_total_sume_layer_25/F");
    outtree->Branch("total_sume_layer_26", &b_total_sume_layer_26, "b_total_sume_layer_26/F");
    outtree->Branch("total_sume_layer_27", &b_total_sume_layer_27, "b_total_sume_layer_27/F");
    outtree->Branch("total_sume_layer_28", &b_total_sume_layer_28, "b_total_sume_layer_28/F");
    outtree->Branch("total_sume_layer_29", &b_total_sume_layer_29, "b_total_sume_layer_29/F");
    outtree->Branch("total_sume_layer_30", &b_total_sume_layer_30, "b_total_sume_layer_30/F");
    outtree->Branch("total_sume_layer_31", &b_total_sume_layer_31, "b_total_sume_layer_31/F");
    outtree->Branch("total_sume_layer_32", &b_total_sume_layer_32, "b_total_sume_layer_32/F");
    outtree->Branch("total_sume_layer_33", &b_total_sume_layer_33, "b_total_sume_layer_33/F");
    outtree->Branch("total_sume_layer_34", &b_total_sume_layer_34, "b_total_sume_layer_34/F");
    outtree->Branch("total_sume_layer_35", &b_total_sume_layer_35, "b_total_sume_layer_35/F");
    outtree->Branch("total_sume_layer_36", &b_total_sume_layer_36, "b_total_sume_layer_36/F");
    outtree->Branch("total_sume_layer_37", &b_total_sume_layer_37, "b_total_sume_layer_37/F");
    outtree->Branch("total_sume_layer_38", &b_total_sume_layer_38, "b_total_sume_layer_38/F");
    outtree->Branch("total_sume_layer_39", &b_total_sume_layer_39, "b_total_sume_layer_39/F");
    outtree->Branch("total_sume_layer_40", &b_total_sume_layer_40, "b_total_sume_layer_40/F");
    outtree->Branch("total_sume_layer_41", &b_total_sume_layer_41, "b_total_sume_layer_41/F");
    outtree->Branch("total_sume_layer_42", &b_total_sume_layer_42, "b_total_sume_layer_42/F");
    outtree->Branch("total_sume_layer_43", &b_total_sume_layer_43, "b_total_sume_layer_43F");
    outtree->Branch("total_sume_layer_44", &b_total_sume_layer_44, "b_total_sume_layer_44F");
    outtree->Branch("total_sume_layer_45", &b_total_sume_layer_45, "b_total_sume_layer_45F");
    outtree->Branch("total_sume_layer_46", &b_total_sume_layer_46, "b_total_sume_layer_46F");
    outtree->Branch("total_sume_layer_47", &b_total_sume_layer_47, "b_total_sume_layer_47F");
    outtree->Branch("total_sume_layer_48", &b_total_sume_layer_48, "b_total_sume_layer_48F");
    outtree->Branch("total_sume_layer_49", &b_total_sume_layer_49, "b_total_sume_layer_49F");
    outtree->Branch("total_sume_layer_50", &b_total_sume_layer_50, "b_total_sume_layer_50/F");
    outtree->Branch("total_sume_layer_51", &b_total_sume_layer_51, "b_total_sume_layer_51/F");
    outtree->Branch("total_sume_layer_52", &b_total_sume_layer_52, "b_total_sume_layer_52/F");
    outtree->Branch("total_sume_layer_53", &b_total_sume_layer_53, "b_total_sume_layer_53/F");
    outtree->Branch("total_sume_layer_54", &b_total_sume_layer_54, "b_total_sume_layer_54/F");
    outtree->Branch("total_sume_layer_55", &b_total_sume_layer_55, "b_total_sume_layer_55/F");
    outtree->Branch("total_sume_layer_n_0", &b_total_sume_layer_n_0, "b_total_sume_layer_n_0/F");
    outtree->Branch("total_sume_layer_n_1", &b_total_sume_layer_n_1, "b_total_sume_layer_n_1/F");
    outtree->Branch("total_sume_layer_n_2", &b_total_sume_layer_n_2, "b_total_sume_layer_n_2/F");
    outtree->Branch("total_sume_layer_n_3", &b_total_sume_layer_n_3, "b_total_sume_layer_n_3/F");
    outtree->Branch("total_sume_layer_n_4", &b_total_sume_layer_n_4, "b_total_sume_layer_n_4/F");
    outtree->Branch("total_sume_layer_n_5", &b_total_sume_layer_n_5, "b_total_sume_layer_n_5/F");
    outtree->Branch("total_sume_layer_n_6", &b_total_sume_layer_n_6, "b_total_sume_layer_n_6/F");
    outtree->Branch("total_sume_layer_n_7", &b_total_sume_layer_n_7, "b_total_sume_layer_n_7/F");
    outtree->Branch("total_sume_layer_n_8", &b_total_sume_layer_n_8, "b_total_sume_layer_n_8/F");
    outtree->Branch("total_sume_layer_n_9", &b_total_sume_layer_n_9, "b_total_sume_layer_n_9/F");
    outtree->Branch("total_sume_layer_n_10", &b_total_sume_layer_n_10, "b_total_sume_layer_n_10/F");
    outtree->Branch("total_sume_layer_n_11", &b_total_sume_layer_n_11, "b_total_sume_layer_n_11/F");
    outtree->Branch("total_sume_layer_n_12", &b_total_sume_layer_n_12, "b_total_sume_layer_n_12/F");
    outtree->Branch("total_sume_layer_n_13", &b_total_sume_layer_n_13, "b_total_sume_layer_n_13/F");
    outtree->Branch("total_sume_layer_n_14", &b_total_sume_layer_n_14, "b_total_sume_layer_n_14/F");
    outtree->Branch("total_sume_layer_n_15", &b_total_sume_layer_n_15, "b_total_sume_layer_n_15/F");
    outtree->Branch("total_sume_layer_n_16", &b_total_sume_layer_n_16, "b_total_sume_layer_n_16/F");
    outtree->Branch("total_sume_layer_n_17", &b_total_sume_layer_n_17, "b_total_sume_layer_n_17/F");
    outtree->Branch("total_sume_layer_n_18", &b_total_sume_layer_n_18, "b_total_sume_layer_n_18/F");
    outtree->Branch("total_sume_layer_n_19", &b_total_sume_layer_n_19, "b_total_sume_layer_n_19/F");
    outtree->Branch("total_sume_layer_n_20", &b_total_sume_layer_n_20, "b_total_sume_layer_n_20/F");
    outtree->Branch("total_sume_layer_n_21", &b_total_sume_layer_n_21, "b_total_sume_layer_n_21/F");
    outtree->Branch("total_sume_layer_n_22", &b_total_sume_layer_n_22, "b_total_sume_layer_n_22/F");
    outtree->Branch("total_sume_layer_n_23", &b_total_sume_layer_n_23, "b_total_sume_layer_n_23/F");
    outtree->Branch("total_sume_layer_n_24", &b_total_sume_layer_n_24, "b_total_sume_layer_n_24/F");
    outtree->Branch("total_sume_layer_n_25", &b_total_sume_layer_n_25, "b_total_sume_layer_n_25/F");
    outtree->Branch("total_sume_layer_n_26", &b_total_sume_layer_n_26, "b_total_sume_layer_n_26/F");
    outtree->Branch("total_sume_layer_n_27", &b_total_sume_layer_n_27, "b_total_sume_layer_n_27/F");
    outtree->Branch("total_sume_layer_n_28", &b_total_sume_layer_n_28, "b_total_sume_layer_n_28/F");
    outtree->Branch("total_sume_layer_n_29", &b_total_sume_layer_n_29, "b_total_sume_layer_n_29/F");
    outtree->Branch("total_sume_layer_n_30", &b_total_sume_layer_n_30, "b_total_sume_layer_n_30/F");
    outtree->Branch("total_sume_layer_n_31", &b_total_sume_layer_n_31, "b_total_sume_layer_n_31/F");
    outtree->Branch("total_sume_layer_n_32", &b_total_sume_layer_n_32, "b_total_sume_layer_n_32/F");
    outtree->Branch("total_sume_layer_n_33", &b_total_sume_layer_n_33, "b_total_sume_layer_n_33/F");
    outtree->Branch("total_sume_layer_n_34", &b_total_sume_layer_n_34, "b_total_sume_layer_n_34/F");
    outtree->Branch("total_sume_layer_n_35", &b_total_sume_layer_n_35, "b_total_sume_layer_n_35/F");
    outtree->Branch("total_sume_layer_n_36", &b_total_sume_layer_n_36, "b_total_sume_layer_n_36/F");
    outtree->Branch("total_sume_layer_n_37", &b_total_sume_layer_n_37, "b_total_sume_layer_n_37/F");
    outtree->Branch("total_sume_layer_n_38", &b_total_sume_layer_n_38, "b_total_sume_layer_n_38/F");
    outtree->Branch("total_sume_layer_n_39", &b_total_sume_layer_n_39, "b_total_sume_layer_n_39/F");
    outtree->Branch("total_sume_layer_n_40", &b_total_sume_layer_n_40, "b_total_sume_layer_n_40/F");
    outtree->Branch("total_sume_layer_n_41", &b_total_sume_layer_n_41, "b_total_sume_layer_n_41/F");
    outtree->Branch("total_sume_layer_n_42", &b_total_sume_layer_n_42, "b_total_sume_layer_n_42/F");
    outtree->Branch("total_sume_layer_n_43", &b_total_sume_layer_n_43, "b_total_sume_layer_n_43F");
    outtree->Branch("total_sume_layer_n_44", &b_total_sume_layer_n_44, "b_total_sume_layer_n_44F");
    outtree->Branch("total_sume_layer_n_45", &b_total_sume_layer_n_45, "b_total_sume_layer_n_45F");
    outtree->Branch("total_sume_layer_n_46", &b_total_sume_layer_n_46, "b_total_sume_layer_n_46F");
    outtree->Branch("total_sume_layer_n_47", &b_total_sume_layer_n_47, "b_total_sume_layer_n_47F");
    outtree->Branch("total_sume_layer_n_48", &b_total_sume_layer_n_48, "b_total_sume_layer_n_48F");
    outtree->Branch("total_sume_layer_n_49", &b_total_sume_layer_n_49, "b_total_sume_layer_n_49F");
    outtree->Branch("total_sume_layer_n_50", &b_total_sume_layer_n_50, "b_total_sume_layer_n_50/F");
    outtree->Branch("total_sume_layer_n_51", &b_total_sume_layer_n_51, "b_total_sume_layer_n_51/F");
    outtree->Branch("total_sume_layer_n_52", &b_total_sume_layer_n_52, "b_total_sume_layer_n_52/F");
    outtree->Branch("total_sume_layer_n_53", &b_total_sume_layer_n_53, "b_total_sume_layer_n_53/F");
    outtree->Branch("total_sume_layer_n_54", &b_total_sume_layer_n_54, "b_total_sume_layer_n_54/F");
    outtree->Branch("total_sume_layer_n_55", &b_total_sume_layer_n_55, "b_total_sume_layer_n_55/F");
    outtree->Branch("total_weighte_layer_0", &b_total_weighte_layer_0, "b_total_weighte_layer_0/F");
    outtree->Branch("total_weighte_layer_1", &b_total_weighte_layer_1, "b_total_weighte_layer_1/F");
    outtree->Branch("total_weighte_layer_2", &b_total_weighte_layer_2, "b_total_weighte_layer_2/F");
    outtree->Branch("total_weighte_layer_3", &b_total_weighte_layer_3, "b_total_weighte_layer_3/F");
    outtree->Branch("total_weighte_layer_4", &b_total_weighte_layer_4, "b_total_weighte_layer_4/F");
    outtree->Branch("total_weighte_layer_5", &b_total_weighte_layer_5, "b_total_weighte_layer_5/F");
    outtree->Branch("total_weighte_layer_6", &b_total_weighte_layer_6, "b_total_weighte_layer_6/F");
    outtree->Branch("total_weighte_layer_7", &b_total_weighte_layer_7, "b_total_weighte_layer_7/F");
    outtree->Branch("total_weighte_layer_8", &b_total_weighte_layer_8, "b_total_weighte_layer_8/F");
    outtree->Branch("total_weighte_layer_9", &b_total_weighte_layer_9, "b_total_weighte_layer_9/F");
    outtree->Branch("total_weighte_layer_10", &b_total_weighte_layer_10, "b_total_weighte_layer_10/F");
    outtree->Branch("total_weighte_layer_11", &b_total_weighte_layer_11, "b_total_weighte_layer_11/F");
    outtree->Branch("total_weighte_layer_12", &b_total_weighte_layer_12, "b_total_weighte_layer_12/F");
    outtree->Branch("total_weighte_layer_13", &b_total_weighte_layer_13, "b_total_weighte_layer_13/F");
    outtree->Branch("total_weighte_layer_14", &b_total_weighte_layer_14, "b_total_weighte_layer_14/F");
    outtree->Branch("total_weighte_layer_15", &b_total_weighte_layer_15, "b_total_weighte_layer_15/F");
    outtree->Branch("total_weighte_layer_16", &b_total_weighte_layer_16, "b_total_weighte_layer_16/F");
    outtree->Branch("total_weighte_layer_17", &b_total_weighte_layer_17, "b_total_weighte_layer_17/F");
    outtree->Branch("total_weighte_layer_18", &b_total_weighte_layer_18, "b_total_weighte_layer_18/F");
    outtree->Branch("total_weighte_layer_19", &b_total_weighte_layer_19, "b_total_weighte_layer_19/F");
    outtree->Branch("total_weighte_layer_20", &b_total_weighte_layer_20, "b_total_weighte_layer_20/F");
    outtree->Branch("total_weighte_layer_21", &b_total_weighte_layer_21, "b_total_weighte_layer_21/F");
    outtree->Branch("total_weighte_layer_22", &b_total_weighte_layer_22, "b_total_weighte_layer_22/F");
    outtree->Branch("total_weighte_layer_23", &b_total_weighte_layer_23, "b_total_weighte_layer_23/F");
    outtree->Branch("total_weighte_layer_24", &b_total_weighte_layer_24, "b_total_weighte_layer_24/F");
    outtree->Branch("total_weighte_layer_25", &b_total_weighte_layer_25, "b_total_weighte_layer_25/F");
    outtree->Branch("total_weighte_layer_26", &b_total_weighte_layer_26, "b_total_weighte_layer_26/F");
    outtree->Branch("total_weighte_layer_27", &b_total_weighte_layer_27, "b_total_weighte_layer_27/F");
    outtree->Branch("total_weighte_layer_28", &b_total_weighte_layer_28, "b_total_weighte_layer_28/F");
    outtree->Branch("total_weighte_layer_29", &b_total_weighte_layer_29, "b_total_weighte_layer_29/F");
    outtree->Branch("total_weighte_layer_30", &b_total_weighte_layer_30, "b_total_weighte_layer_30/F");
    outtree->Branch("total_weighte_layer_31", &b_total_weighte_layer_31, "b_total_weighte_layer_31/F");
    outtree->Branch("total_weighte_layer_32", &b_total_weighte_layer_32, "b_total_weighte_layer_32/F");
    outtree->Branch("total_weighte_layer_33", &b_total_weighte_layer_33, "b_total_weighte_layer_33/F");
    outtree->Branch("total_weighte_layer_34", &b_total_weighte_layer_34, "b_total_weighte_layer_34/F");
    outtree->Branch("total_weighte_layer_35", &b_total_weighte_layer_35, "b_total_weighte_layer_35/F");
    outtree->Branch("total_weighte_layer_36", &b_total_weighte_layer_36, "b_total_weighte_layer_36/F");
    outtree->Branch("total_weighte_layer_37", &b_total_weighte_layer_37, "b_total_weighte_layer_37/F");
    outtree->Branch("total_weighte_layer_38", &b_total_weighte_layer_38, "b_total_weighte_layer_38/F");
    outtree->Branch("total_weighte_layer_39", &b_total_weighte_layer_39, "b_total_weighte_layer_39/F");
    outtree->Branch("total_weighte_layer_40", &b_total_weighte_layer_40, "b_total_weighte_layer_40/F");
    outtree->Branch("total_weighte_layer_41", &b_total_weighte_layer_41, "b_total_weighte_layer_41/F");
    outtree->Branch("total_weighte_layer_42", &b_total_weighte_layer_42, "b_total_weighte_layer_42/F");
    outtree->Branch("total_weighte_layer_43", &b_total_weighte_layer_43, "b_total_weighte_layer_43F");
    outtree->Branch("total_weighte_layer_44", &b_total_weighte_layer_44, "b_total_weighte_layer_44F");
    outtree->Branch("total_weighte_layer_45", &b_total_weighte_layer_45, "b_total_weighte_layer_45F");
    outtree->Branch("total_weighte_layer_46", &b_total_weighte_layer_46, "b_total_weighte_layer_46F");
    outtree->Branch("total_weighte_layer_47", &b_total_weighte_layer_47, "b_total_weighte_layer_47F");
    outtree->Branch("total_weighte_layer_48", &b_total_weighte_layer_48, "b_total_weighte_layer_48F");
    outtree->Branch("total_weighte_layer_49", &b_total_weighte_layer_49, "b_total_weighte_layer_49F");
    outtree->Branch("total_weighte_layer_50", &b_total_weighte_layer_50, "b_total_weighte_layer_50/F");
    outtree->Branch("total_weighte_layer_51", &b_total_weighte_layer_51, "b_total_weighte_layer_51/F");
    outtree->Branch("total_weighte_layer_52", &b_total_weighte_layer_52, "b_total_weighte_layer_52/F");
    outtree->Branch("total_weighte_layer_53", &b_total_weighte_layer_53, "b_total_weighte_layer_53/F");
    outtree->Branch("total_weighte_layer_54", &b_total_weighte_layer_54, "b_total_weighte_layer_54/F");
    outtree->Branch("total_weighte_layer_55", &b_total_weighte_layer_55, "b_total_weighte_layer_55/F");
    outtree->Branch("total_weighte_layer_n_0", &b_total_weighte_layer_n_0, "b_total_weighte_layer_n_0/F");
    outtree->Branch("total_weighte_layer_n_1", &b_total_weighte_layer_n_1, "b_total_weighte_layer_n_1/F");
    outtree->Branch("total_weighte_layer_n_2", &b_total_weighte_layer_n_2, "b_total_weighte_layer_n_2/F");
    outtree->Branch("total_weighte_layer_n_3", &b_total_weighte_layer_n_3, "b_total_weighte_layer_n_3/F");
    outtree->Branch("total_weighte_layer_n_4", &b_total_weighte_layer_n_4, "b_total_weighte_layer_n_4/F");
    outtree->Branch("total_weighte_layer_n_5", &b_total_weighte_layer_n_5, "b_total_weighte_layer_n_5/F");
    outtree->Branch("total_weighte_layer_n_6", &b_total_weighte_layer_n_6, "b_total_weighte_layer_n_6/F");
    outtree->Branch("total_weighte_layer_n_7", &b_total_weighte_layer_n_7, "b_total_weighte_layer_n_7/F");
    outtree->Branch("total_weighte_layer_n_8", &b_total_weighte_layer_n_8, "b_total_weighte_layer_n_8/F");
    outtree->Branch("total_weighte_layer_n_9", &b_total_weighte_layer_n_9, "b_total_weighte_layer_n_9/F");
    outtree->Branch("total_weighte_layer_n_10", &b_total_weighte_layer_n_10, "b_total_weighte_layer_n_10/F");
    outtree->Branch("total_weighte_layer_n_11", &b_total_weighte_layer_n_11, "b_total_weighte_layer_n_11/F");
    outtree->Branch("total_weighte_layer_n_12", &b_total_weighte_layer_n_12, "b_total_weighte_layer_n_12/F");
    outtree->Branch("total_weighte_layer_n_13", &b_total_weighte_layer_n_13, "b_total_weighte_layer_n_13/F");
    outtree->Branch("total_weighte_layer_n_14", &b_total_weighte_layer_n_14, "b_total_weighte_layer_n_14/F");
    outtree->Branch("total_weighte_layer_n_15", &b_total_weighte_layer_n_15, "b_total_weighte_layer_n_15/F");
    outtree->Branch("total_weighte_layer_n_16", &b_total_weighte_layer_n_16, "b_total_weighte_layer_n_16/F");
    outtree->Branch("total_weighte_layer_n_17", &b_total_weighte_layer_n_17, "b_total_weighte_layer_n_17/F");
    outtree->Branch("total_weighte_layer_n_18", &b_total_weighte_layer_n_18, "b_total_weighte_layer_n_18/F");
    outtree->Branch("total_weighte_layer_n_19", &b_total_weighte_layer_n_19, "b_total_weighte_layer_n_19/F");
    outtree->Branch("total_weighte_layer_n_20", &b_total_weighte_layer_n_20, "b_total_weighte_layer_n_20/F");
    outtree->Branch("total_weighte_layer_n_21", &b_total_weighte_layer_n_21, "b_total_weighte_layer_n_21/F");
    outtree->Branch("total_weighte_layer_n_22", &b_total_weighte_layer_n_22, "b_total_weighte_layer_n_22/F");
    outtree->Branch("total_weighte_layer_n_23", &b_total_weighte_layer_n_23, "b_total_weighte_layer_n_23/F");
    outtree->Branch("total_weighte_layer_n_24", &b_total_weighte_layer_n_24, "b_total_weighte_layer_n_24/F");
    outtree->Branch("total_weighte_layer_n_25", &b_total_weighte_layer_n_25, "b_total_weighte_layer_n_25/F");
    outtree->Branch("total_weighte_layer_n_26", &b_total_weighte_layer_n_26, "b_total_weighte_layer_n_26/F");
    outtree->Branch("total_weighte_layer_n_27", &b_total_weighte_layer_n_27, "b_total_weighte_layer_n_27/F");
    outtree->Branch("total_weighte_layer_n_28", &b_total_weighte_layer_n_28, "b_total_weighte_layer_n_28/F");
    outtree->Branch("total_weighte_layer_n_29", &b_total_weighte_layer_n_29, "b_total_weighte_layer_n_29/F");
    outtree->Branch("total_weighte_layer_n_30", &b_total_weighte_layer_n_30, "b_total_weighte_layer_n_30/F");
    outtree->Branch("total_weighte_layer_n_31", &b_total_weighte_layer_n_31, "b_total_weighte_layer_n_31/F");
    outtree->Branch("total_weighte_layer_n_32", &b_total_weighte_layer_n_32, "b_total_weighte_layer_n_32/F");
    outtree->Branch("total_weighte_layer_n_33", &b_total_weighte_layer_n_33, "b_total_weighte_layer_n_33/F");
    outtree->Branch("total_weighte_layer_n_34", &b_total_weighte_layer_n_34, "b_total_weighte_layer_n_34/F");
    outtree->Branch("total_weighte_layer_n_35", &b_total_weighte_layer_n_35, "b_total_weighte_layer_n_35/F");
    outtree->Branch("total_weighte_layer_n_36", &b_total_weighte_layer_n_36, "b_total_weighte_layer_n_36/F");
    outtree->Branch("total_weighte_layer_n_37", &b_total_weighte_layer_n_37, "b_total_weighte_layer_n_37/F");
    outtree->Branch("total_weighte_layer_n_38", &b_total_weighte_layer_n_38, "b_total_weighte_layer_n_38/F");
    outtree->Branch("total_weighte_layer_n_39", &b_total_weighte_layer_n_39, "b_total_weighte_layer_n_39/F");
    outtree->Branch("total_weighte_layer_n_40", &b_total_weighte_layer_n_40, "b_total_weighte_layer_n_40/F");
    outtree->Branch("total_weighte_layer_n_41", &b_total_weighte_layer_n_41, "b_total_weighte_layer_n_41/F");
    outtree->Branch("total_weighte_layer_n_42", &b_total_weighte_layer_n_42, "b_total_weighte_layer_n_42/F");
    outtree->Branch("total_weighte_layer_n_43", &b_total_weighte_layer_n_43, "b_total_weighte_layer_n_43F");
    outtree->Branch("total_weighte_layer_n_44", &b_total_weighte_layer_n_44, "b_total_weighte_layer_n_44F");
    outtree->Branch("total_weighte_layer_n_45", &b_total_weighte_layer_n_45, "b_total_weighte_layer_n_45F");
    outtree->Branch("total_weighte_layer_n_46", &b_total_weighte_layer_n_46, "b_total_weighte_layer_n_46F");
    outtree->Branch("total_weighte_layer_n_47", &b_total_weighte_layer_n_47, "b_total_weighte_layer_n_47F");
    outtree->Branch("total_weighte_layer_n_48", &b_total_weighte_layer_n_48, "b_total_weighte_layer_n_48F");
    outtree->Branch("total_weighte_layer_n_49", &b_total_weighte_layer_n_49, "b_total_weighte_layer_n_49F");
    outtree->Branch("total_weighte_layer_n_50", &b_total_weighte_layer_n_50, "b_total_weighte_layer_n_50/F");
    outtree->Branch("total_weighte_layer_n_51", &b_total_weighte_layer_n_51, "b_total_weighte_layer_n_51/F");
    outtree->Branch("total_weighte_layer_n_52", &b_total_weighte_layer_n_52, "b_total_weighte_layer_n_52/F");
    outtree->Branch("total_weighte_layer_n_53", &b_total_weighte_layer_n_53, "b_total_weighte_layer_n_53/F");
    outtree->Branch("total_weighte_layer_n_54", &b_total_weighte_layer_n_54, "b_total_weighte_layer_n_54/F");
    outtree->Branch("total_weighte_layer_n_55", &b_total_weighte_layer_n_55, "b_total_weighte_layer_n_55/F");

    outtree->Branch("nhit_ratio", &b_nhit_ratio , "b_nhit_ratio/F");
    outtree->Branch("sume_ratio", &b_sume_ratio , "b_sume_ratio/F");
    outtree->Branch("weighte_ratio", &b_weighte_ratio , "b_weighte_ratio/F");

    // Begin of reading files and writing trees
    for (int i_energy = 0; i_energy < N_ENERGIES; i_energy++) {
        cout << "Energy = " << energies[i_energy] << " GeV" << endl;
        // Read files for the 3 parts of the detector
	string filename_ecal = filenames_ecal[i_energy];
        TFile *file_ecal = new TFile(filename_ecal.c_str(), "read");
        TTree *tree_ecal = (TTree*)file_ecal->Get("ecal");
	string filename_hcal = filenames_hcal[i_energy];
        TFile *file_hcal = new TFile(filename_hcal.c_str(), "read");
        TTree *tree_hcal = (TTree*)file_hcal->Get("hcal");

	tree_ecal->AddFriend(tree_hcal);

        Long64_t nentries = tree_ecal->GetEntries();

	cout << "nentries (ecal, ahcal): (" <<tree_ecal->GetEntries()<<", "<<tree_hcal->GetEntries()<<")"<<endl;
	
        vector<float> *ecal_hit_energy = 0;
        vector<float> *ecal_hit_x = 0;
        vector<float> *ecal_hit_y = 0;
        vector<float> *ecal_hit_z = 0;
        vector<int> *ecal_hit_isMasked = 0;
        vector<int> *ecal_hit_slab = 0;
        
        TBranch *ecal_bhit_energy = 0;
        TBranch *ecal_bhit_x = 0;
        TBranch *ecal_bhit_y = 0;
        TBranch *ecal_bhit_z = 0;
        TBranch *ecal_bhit_isMasked = 0;
        TBranch *ecal_bhit_slab= 0;

        tree_ecal->SetBranchAddress("hit_energy", &ecal_hit_energy, &ecal_bhit_energy);
        tree_ecal->SetBranchAddress("hit_x", &ecal_hit_x, &ecal_bhit_x);
        tree_ecal->SetBranchAddress("hit_y", &ecal_hit_y, &ecal_bhit_y);
        tree_ecal->SetBranchAddress("hit_z", &ecal_hit_z, &ecal_bhit_z);
        tree_ecal->SetBranchAddress("hit_isMasked", &ecal_hit_isMasked, &ecal_bhit_isMasked);
        tree_ecal->SetBranchAddress("hit_slab", &ecal_hit_slab, &ecal_bhit_slab);
        	
	vector<float> *hcal_hit_energy = 0;
        vector<float> *hcal_hit_x = 0;
        vector<float> *hcal_hit_y = 0;
        vector<float> *hcal_hit_z = 0;
        vector<int> *hcal_hit_isMasked = 0;
        vector<int> *hcal_hit_slab = 0;

        TBranch *hcal_bhit_energy = 0;
        TBranch *hcal_bhit_x = 0;
        TBranch *hcal_bhit_y = 0;
        TBranch *hcal_bhit_z = 0;
        TBranch *hcal_bhit_isMasked = 0;
        TBranch *hcal_bhit_slab= 0;

        tree_hcal->SetBranchAddress("hit_energy", &hcal_hit_energy, &hcal_bhit_energy);
        tree_hcal->SetBranchAddress("hit_x", &hcal_hit_x, &hcal_bhit_x);
        tree_hcal->SetBranchAddress("hit_y", &hcal_hit_y, &hcal_bhit_y);
        tree_hcal->SetBranchAddress("hit_z", &hcal_hit_z, &hcal_bhit_z);
        tree_hcal->SetBranchAddress("hit_isMasked", &hcal_hit_isMasked, &hcal_bhit_isMasked);
        tree_hcal->SetBranchAddress("hit_slab", &hcal_hit_slab, &hcal_bhit_slab);

	for (int i_event = 0; i_event < (int)nentries; i_event++) {
	//for (int i_event = 0; i_event < 100; i_event++) {
	  // Declarations
	  int ecal_nhit = 0;     
	  float ecal_sume = 0.;   
	  float ecal_weighte = 0.; 
	  float ecal_mol_value = 0.;
	  float ecal_bar_xyzr[4];
	  float ecal_nhit_layer_array[N_ECAL_LAYERS];
	  float ecal_nhit_layer_n_array[N_ECAL_LAYERS];
	  float ecal_sume_layer_array[N_ECAL_LAYERS];
          float ecal_sume_layer_n_array[N_ECAL_LAYERS];
	  float ecal_weighte_layer_array[N_ECAL_LAYERS];
          float ecal_weighte_layer_n_array[N_ECAL_LAYERS];
	  float ecal_bar_layer_array[N_ECAL_LAYERS][3];
	  float ecal_radius90_layer_array[N_ECAL_LAYERS];

	  int hcal_nhit = 0;
          float hcal_sume = 0.;
          float hcal_weighte = 0.;
          float hcal_mol_value = 0.;
          float hcal_bar_xyzr[4];
          float hcal_nhit_layer_array[N_HCAL_LAYERS];
          float hcal_nhit_layer_n_array[N_HCAL_LAYERS];
          float hcal_sume_layer_array[N_HCAL_LAYERS];
          float hcal_sume_layer_n_array[N_HCAL_LAYERS];
          float hcal_weighte_layer_array[N_HCAL_LAYERS];
          float hcal_weighte_layer_n_array[N_HCAL_LAYERS];
          float hcal_bar_layer_array[N_HCAL_LAYERS][3];
          float hcal_radius90_layer_array[N_HCAL_LAYERS];

	  int total_nhit = 0;
          float total_sume = 0;
          float total_weighte = 0;
          float total_mol_value = 0.;
          float total_bar_xyzr[4];
          float total_nhit_layer_array[N_TOTAL_LAYERS];
          float total_nhit_layer_n_array[N_TOTAL_LAYERS];
          float total_sume_layer_array[N_TOTAL_LAYERS];
          float total_sume_layer_n_array[N_TOTAL_LAYERS];
          float total_weighte_layer_array[N_TOTAL_LAYERS];
          float total_weighte_layer_n_array[N_TOTAL_LAYERS];
          float total_bar_layer_array[N_TOTAL_LAYERS][3];
          float total_radius90_layer_array[N_TOTAL_LAYERS];

	  for(int i = 0; i<N_ECAL_LAYERS;i++) {
            ecal_nhit_layer_array[i] = 0.;
            ecal_nhit_layer_n_array[i] = 0.;
            ecal_sume_layer_array[i] = 0.;
            ecal_sume_layer_n_array[i] = 0.;
            ecal_weighte_layer_array[i] = 0.;
            ecal_weighte_layer_n_array[i] = 0.;
            ecal_bar_layer_array[i][0] = 0.;
            ecal_bar_layer_array[i][1] = 0.;
            ecal_bar_layer_array[i][2] = 0.;
            ecal_radius90_layer_array[i] = 0.;
	  }
	  for(int i = 0; i<N_HCAL_LAYERS;i++) {
	    hcal_nhit_layer_array[i] = 0.;
            hcal_nhit_layer_n_array[i] = 0.;
            hcal_sume_layer_array[i] = 0.;
            hcal_sume_layer_n_array[i] = 0.;
            hcal_weighte_layer_array[i] = 0.;
            hcal_weighte_layer_n_array[i] = 0.;
            hcal_bar_layer_array[i][0] = 0.;
            hcal_bar_layer_array[i][1] = 0.;
            hcal_bar_layer_array[i][2] = 0.;
            hcal_radius90_layer_array[i] = 0.;
	  }
	  for(int i = 0; i<N_TOTAL_LAYERS;i++) {
	    total_nhit_layer_array[i] = 0.;
            total_nhit_layer_n_array[i] = 0.;
            total_sume_layer_array[i] = 0.;
            total_sume_layer_n_array[i] = 0.;
            total_weighte_layer_array[i] = 0.;
            total_weighte_layer_n_array[i] = 0.;
            total_bar_layer_array[i][0] = 0.;
            total_bar_layer_array[i][1] = 0.;
            total_bar_layer_array[i][2] = 0.;
            total_radius90_layer_array[i] = 0.;
	  }

	  // Get values
	  tree_ecal->GetEntry(i_event);
	  tree_hcal->GetEntry(i_event);
	  
	  if(i_event % 1000 == 0) cout << "Event " << to_string(i_event) << endl;
	  
	  get_energy_ecal(ecal_nhit, ecal_sume, ecal_weighte, ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_isMasked, masked);
	  barycenter_ecal(ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_x, ecal_hit_y, ecal_hit_z, ecal_bar_xyzr, ecal_hit_isMasked, masked);
	  hits_layer_ecal(ecal_nhit_layer_array, ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_isMasked, masked, false, "nhit");
	  hits_layer_ecal(ecal_nhit_layer_n_array, ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_isMasked, masked, true, "nhit");
          hits_layer_ecal(ecal_sume_layer_array, ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_isMasked, masked, false, "sume");
          hits_layer_ecal(ecal_sume_layer_n_array, ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_isMasked, masked, true, "sume");
	  hits_layer_ecal(ecal_weighte_layer_array, ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_isMasked, masked, false, "weight");
          hits_layer_ecal(ecal_weighte_layer_n_array, ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_isMasked, masked, true, "weight");
	  
	  float ecal_hits_max_distance_value = hits_max_distance(ecal_hit_slab, ecal_hit_x, ecal_hit_x, ecal_hit_isMasked, masked);
	  float ecal_MIP_Likeness_value = MIP_Likeness_ecal(ecal_nhit_layer_array);
	  bool ecal_shower_bool = is_Shower_ecal(ecal_sume, ecal_sume_layer_array);

	  float ecal_nhit_shower_maxvalue = 0.;
	  float ecal_nhit_shower_maxvalue_n = 0.;
	  int ecal_nhit_ilayermax = -1;
	  int ecal_nhit_ilayerstart = -1;
	  int ecal_nhit_ilayerend = -1;
	  int ecal_nhit_ilayerstart_10 = -1;
          int ecal_nhit_ilayerend_10 = -1;
	  float ecal_nhit_shower_averagevalue = ecal_nhit/N_ECAL_LAYERS;
	  shower_variables_ecal(ecal_nhit, ecal_nhit_layer_array, ecal_nhit_layer_n_array, ecal_nhit_shower_maxvalue, ecal_nhit_shower_maxvalue_n, ecal_nhit_ilayermax,
			   ecal_nhit_ilayerstart, ecal_nhit_ilayerstart_10, ecal_nhit_ilayerend, ecal_nhit_ilayerend_10, "ecal_nhit", ecal_shower_bool);

          float ecal_sume_shower_maxvalue = 0.;
          float ecal_sume_shower_maxvalue_n = 0.;
          int ecal_sume_ilayermax = -1;
          int ecal_sume_ilayerstart = -1;
          int ecal_sume_ilayerend = -1;
          int ecal_sume_ilayerstart_10 = -1;
          int ecal_sume_ilayerend_10 = -1;
          float ecal_sume_shower_averagevalue = ecal_sume/N_ECAL_LAYERS;
          shower_variables_ecal(ecal_sume, ecal_sume_layer_array, ecal_sume_layer_n_array, ecal_sume_shower_maxvalue, ecal_sume_shower_maxvalue_n, ecal_sume_ilayermax,
                           ecal_sume_ilayerstart, ecal_sume_ilayerstart_10, ecal_sume_ilayerend, ecal_sume_ilayerend_10, "ecal_sume", ecal_shower_bool);
         
          float ecal_weighte_shower_maxvalue = 0.;
	  float ecal_weighte_shower_maxvalue_n = 0.;
          int ecal_weighte_ilayermax = -1;
          int ecal_weighte_ilayerstart = -1;
          int ecal_weighte_ilayerend = -1;
          int ecal_weighte_ilayerstart_10 = -1;
          int ecal_weighte_ilayerend_10 = -1;
	  float ecal_weighte_shower_averagevalue = ecal_weighte/N_ECAL_LAYERS;
	  shower_variables_ecal(ecal_weighte, ecal_weighte_layer_array, ecal_weighte_layer_n_array, ecal_weighte_shower_maxvalue, ecal_weighte_shower_maxvalue_n, ecal_weighte_ilayermax,
                           ecal_weighte_ilayerstart, ecal_weighte_ilayerstart_10, ecal_weighte_ilayerend, ecal_weighte_ilayerend_10, "ecal_weight", ecal_shower_bool);

          ecal_mol_value = moliere_ecal(ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_x, ecal_hit_y, ecal_hit_z, ecal_hit_isMasked, masked, 0.9, ecal_shower_bool);
          radius_layer_ecal(ecal_radius90_layer_array, ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_x, ecal_hit_y, ecal_hit_z, ecal_hit_isMasked, masked, 0.9, ecal_shower_bool);
	  bary_layer_ecal(ecal_bar_layer_array, ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_x, ecal_hit_y, ecal_hit_z, ecal_hit_isMasked, masked);



	  get_energy_hcal(hcal_nhit, hcal_sume, hcal_weighte, hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_isMasked, masked);
          barycenter_hcal(hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_x, hcal_hit_y, hcal_hit_z, hcal_bar_xyzr, hcal_hit_isMasked, masked);
          hits_layer_hcal(hcal_nhit_layer_array, hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_isMasked, masked, false, "nhit");
          hits_layer_hcal(hcal_nhit_layer_n_array, hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_isMasked, masked, true, "nhit");
          hits_layer_hcal(hcal_sume_layer_array, hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_isMasked, masked, false, "sume");
          hits_layer_hcal(hcal_sume_layer_n_array, hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_isMasked, masked, true, "sume");
          hits_layer_hcal(hcal_weighte_layer_array, hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_isMasked, masked, false, "weight");
          hits_layer_hcal(hcal_weighte_layer_n_array, hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_isMasked, masked, true, "weight");

	  //float hcal_hits_max_distance_value = 0.;
	  float hcal_hits_max_distance_value = hits_max_distance(hcal_hit_slab, hcal_hit_x, hcal_hit_x, hcal_hit_isMasked, masked);
          float hcal_MIP_Likeness_value = MIP_Likeness_hcal(hcal_nhit_layer_array);
          bool hcal_shower_bool = is_Shower_hcal(hcal_sume, hcal_sume_layer_array);

          float hcal_nhit_shower_maxvalue = 0.;
          float hcal_nhit_shower_maxvalue_n = 0.;
          int hcal_nhit_ilayermax = -1;
          int hcal_nhit_ilayerstart = -1;
          int hcal_nhit_ilayerend = -1;
          int hcal_nhit_ilayerstart_10 = -1;
          int hcal_nhit_ilayerend_10 = -1;
          float hcal_nhit_shower_averagevalue = hcal_nhit/N_HCAL_LAYERS;
          shower_variables_hcal(hcal_nhit, hcal_nhit_layer_array, hcal_nhit_layer_n_array, hcal_nhit_shower_maxvalue, hcal_nhit_shower_maxvalue_n, hcal_nhit_ilayermax,
				hcal_nhit_ilayerstart, hcal_nhit_ilayerstart_10, hcal_nhit_ilayerend, hcal_nhit_ilayerend_10, "nhit", hcal_shower_bool);

          float hcal_sume_shower_maxvalue = 0.;
          float hcal_sume_shower_maxvalue_n = 0.;
          int hcal_sume_ilayermax = -1;
          int hcal_sume_ilayerstart = -1;
          int hcal_sume_ilayerend = -1;
          int hcal_sume_ilayerstart_10 = -1;
          int hcal_sume_ilayerend_10 = -1;
          float hcal_sume_shower_averagevalue = hcal_sume/N_HCAL_LAYERS;
          shower_variables_hcal(hcal_sume, hcal_sume_layer_array, hcal_sume_layer_n_array, hcal_sume_shower_maxvalue, hcal_sume_shower_maxvalue_n, hcal_sume_ilayermax,
				hcal_sume_ilayerstart, hcal_sume_ilayerstart_10, hcal_sume_ilayerend, hcal_sume_ilayerend_10, "sume", hcal_shower_bool);
          float hcal_weighte_shower_maxvalue = 0.;
          float hcal_weighte_shower_maxvalue_n = 0.;
          int hcal_weighte_ilayermax = -1;
          int hcal_weighte_ilayerstart = -1;
          int hcal_weighte_ilayerend = -1;
          int hcal_weighte_ilayerstart_10 = -1;
          int hcal_weighte_ilayerend_10 = -1;
          float hcal_weighte_shower_averagevalue = hcal_weighte/N_HCAL_LAYERS;
          shower_variables_hcal(hcal_weighte, hcal_weighte_layer_array, hcal_weighte_layer_n_array, hcal_weighte_shower_maxvalue, hcal_weighte_shower_maxvalue_n, hcal_weighte_ilayermax,
				hcal_weighte_ilayerstart, hcal_weighte_ilayerstart_10, hcal_weighte_ilayerend, hcal_weighte_ilayerend_10, "weight", hcal_shower_bool);

          hcal_mol_value = moliere_hcal(hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_x, hcal_hit_y, hcal_hit_z, hcal_hit_isMasked, masked, 0.9, hcal_shower_bool);
          radius_layer_hcal(hcal_radius90_layer_array, hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_x, hcal_hit_y, hcal_hit_z, hcal_hit_isMasked, masked, 0.9, hcal_shower_bool);
          bary_layer_hcal(hcal_bar_layer_array, hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_x, hcal_hit_y, hcal_hit_z, hcal_hit_isMasked, masked);


	  //cout<<"debug 2 trees: "<<ecal_nhit<<" "<<hcal_nhit<<endl;

	  bool total_shower_bool = 0;
	  if((ecal_shower_bool==1) or (hcal_shower_bool==1)) total_shower_bool=1;
	    
	  total_mol_value = moliere_total(ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_x, ecal_hit_y, ecal_hit_z, ecal_hit_isMasked, 
					  hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_x, hcal_hit_y, hcal_hit_z, hcal_hit_isMasked,
					  masked, 0.9, total_shower_bool);

          barycenter_total(ecal_hit_energy, ecal_hit_slab, W_thicknesses, ecal_hit_x, ecal_hit_y, ecal_hit_z, ecal_hit_isMasked,
			   hcal_hit_energy, hcal_hit_slab, S_thicknesses, hcal_hit_x, hcal_hit_y, hcal_hit_z, hcal_hit_isMasked,
			   total_bar_xyzr, masked);
	  
	  total_nhit = ecal_nhit + hcal_nhit;
	  total_sume = ecal_sume + hcal_sume;
	  total_weighte = ecal_weighte + hcal_weighte;
	  
	  for(int i = 0; i<N_TOTAL_LAYERS;i++) {
            if(i<N_ECAL_LAYERS) {
              total_nhit_layer_array[i] = ecal_nhit_layer_array[i];
              total_nhit_layer_n_array[i] = ecal_nhit_layer_n_array[i];
              total_sume_layer_array[i] = ecal_sume_layer_array[i];
              total_sume_layer_n_array[i] = ecal_sume_layer_n_array[i];
              total_weighte_layer_array[i] = ecal_weighte_layer_array[i];
              total_weighte_layer_n_array[i] = ecal_weighte_layer_n_array[i];
              total_bar_layer_array[i][0] = ecal_bar_layer_array[i][0];
              total_bar_layer_array[i][1] = ecal_bar_layer_array[i][1];
              total_bar_layer_array[i][2] = ecal_bar_layer_array[i][2];
              total_radius90_layer_array[i] = ecal_radius90_layer_array[i];
            }
            else{
              total_nhit_layer_array[i] = hcal_nhit_layer_array[i-N_ECAL_LAYERS];
              total_nhit_layer_n_array[i] = hcal_nhit_layer_n_array[i-N_ECAL_LAYERS];
              total_sume_layer_array[i] = hcal_sume_layer_array[i-N_ECAL_LAYERS];
              total_sume_layer_n_array[i] = hcal_sume_layer_n_array[i-N_ECAL_LAYERS];
              total_weighte_layer_array[i] = hcal_weighte_layer_array[i-N_ECAL_LAYERS];
              total_weighte_layer_n_array[i] = hcal_weighte_layer_n_array[i-N_ECAL_LAYERS];
              total_bar_layer_array[i][0] = hcal_bar_layer_array[i-N_ECAL_LAYERS][0];
              total_bar_layer_array[i][1] = hcal_bar_layer_array[i-N_ECAL_LAYERS][1];
              total_bar_layer_array[i][2] = hcal_bar_layer_array[i-N_ECAL_LAYERS][2];
              total_radius90_layer_array[i] = hcal_radius90_layer_array[i-N_ECAL_LAYERS];
            }
          }
	  
	  
	  float total_hits_max_distance_value = 0.;
	  if(ecal_hits_max_distance_value > hcal_hits_max_distance_value ) total_hits_max_distance_value = ecal_hits_max_distance_value;
	  if(ecal_hits_max_distance_value < hcal_hits_max_distance_value ) total_hits_max_distance_value = hcal_hits_max_distance_value;
	  float total_MIP_Likeness_value = MIP_Likeness_total(ecal_nhit_layer_array, hcal_nhit_layer_array);
	  float total_nhit_shower_maxvalue = 0.;
          float total_nhit_shower_maxvalue_n = 0.;
          int total_nhit_ilayermax = -1;
          int total_nhit_ilayerstart = -1;
          int total_nhit_ilayerend = -1;
          int total_nhit_ilayerstart_10 = -1;
          int total_nhit_ilayerend_10 = -1;
          float total_nhit_shower_averagevalue = (ecal_nhit+hcal_nhit)/N_TOTAL_LAYERS;
	  shower_variables_total(total_nhit, total_nhit_layer_array, total_nhit_layer_n_array, total_nhit_shower_maxvalue, total_nhit_shower_maxvalue_n, total_nhit_ilayermax,
                                total_nhit_ilayerstart, total_nhit_ilayerstart_10, total_nhit_ilayerend, total_nhit_ilayerend_10, "nhit", total_shower_bool);

	  float total_sume_shower_maxvalue = 0.;
          float total_sume_shower_maxvalue_n = 0.;
          int total_sume_ilayermax = -1;
          int total_sume_ilayerstart = -1;
          int total_sume_ilayerend = -1;
          int total_sume_ilayerstart_10 = -1;
          int total_sume_ilayerend_10 = -1;
          float total_sume_shower_averagevalue = (ecal_sume+hcal_sume)/N_TOTAL_LAYERS;
	  shower_variables_total(total_sume, total_sume_layer_array, total_sume_layer_n_array, total_sume_shower_maxvalue, total_sume_shower_maxvalue_n, total_sume_ilayermax,
				 total_sume_ilayerstart, total_sume_ilayerstart_10, total_sume_ilayerend, total_sume_ilayerend_10, "sume", total_shower_bool);

	  float total_weighte_shower_maxvalue = 0.;
          float total_weighte_shower_maxvalue_n = 0.;
          int total_weighte_ilayermax = -1;
          int total_weighte_ilayerstart = -1;
          int total_weighte_ilayerend = -1;
          int total_weighte_ilayerstart_10 = -1;
          int total_weighte_ilayerend_10 = -1;
          float total_weighte_shower_averagevalue = (ecal_weighte+hcal_weighte)/N_TOTAL_LAYERS;
	  shower_variables_total(total_weighte, total_weighte_layer_array, total_weighte_layer_n_array, total_weighte_shower_maxvalue, total_weighte_shower_maxvalue_n, total_weighte_ilayermax,
				 total_weighte_ilayerstart, total_weighte_ilayerstart_10, total_weighte_ilayerend, total_weighte_ilayerend_10, "weighte", total_shower_bool);

	  float ratio_nhit = 0.;
	  float ratio_sume = 0.;
	  float ratio_weighte = 0.;
	  
	  hcal_ecal_ratio(ratio_nhit, ratio_sume, ratio_weighte, ecal_nhit, ecal_sume, ecal_weighte, hcal_nhit, hcal_sume, hcal_weighte);

	  is_interaction(b_ecal_interaction, b_hcal_interaction, b_total_interaction, ecal_nhit, hcal_nhit);

	  // Filling the tree
	  //cout<<"ecal nhit: "<<ecal_nhit<<endl;
	  b_ecal_nhit = ecal_nhit;
	  b_ecal_sume = ecal_sume;
	  b_ecal_weighte = ecal_weighte;
	  b_ecal_bar_x = ecal_bar_xyzr[0];
	  b_ecal_bar_y = ecal_bar_xyzr[1];
	  b_ecal_bar_z = ecal_bar_xyzr[2];
	  b_ecal_bar_r = ecal_bar_xyzr[3];
	  b_ecal_mol = ecal_mol_value;
	  b_ecal_hits_max_distance = ecal_hits_max_distance_value;
	  b_ecal_MIP_Likeness = ecal_MIP_Likeness_value;

	  b_ecal_radius90_layer_0 = ecal_radius90_layer_array[0];
          b_ecal_radius90_layer_1 = ecal_radius90_layer_array[1];
          b_ecal_radius90_layer_2 = ecal_radius90_layer_array[2];
          b_ecal_radius90_layer_3 = ecal_radius90_layer_array[3];
          b_ecal_radius90_layer_4 = ecal_radius90_layer_array[4];
          b_ecal_radius90_layer_5 = ecal_radius90_layer_array[5];
          b_ecal_radius90_layer_6 = ecal_radius90_layer_array[6];
          b_ecal_radius90_layer_7 = ecal_radius90_layer_array[7];
          b_ecal_radius90_layer_8 = ecal_radius90_layer_array[8];
          b_ecal_radius90_layer_9 = ecal_radius90_layer_array[9];
          b_ecal_radius90_layer_10 = ecal_radius90_layer_array[10];
          b_ecal_radius90_layer_11 = ecal_radius90_layer_array[11];
          b_ecal_radius90_layer_12 = ecal_radius90_layer_array[12];
          b_ecal_radius90_layer_13 = ecal_radius90_layer_array[13];
          b_ecal_radius90_layer_14 = ecal_radius90_layer_array[14];
	  
	  b_ecal_nhit_layer_0 = ecal_nhit_layer_array[0];
	  b_ecal_nhit_layer_1 = ecal_nhit_layer_array[1];
	  b_ecal_nhit_layer_2 = ecal_nhit_layer_array[2];
	  b_ecal_nhit_layer_3 = ecal_nhit_layer_array[3];
	  b_ecal_nhit_layer_4 = ecal_nhit_layer_array[4];
	  b_ecal_nhit_layer_5 = ecal_nhit_layer_array[5];
	  b_ecal_nhit_layer_6 = ecal_nhit_layer_array[6];
	  b_ecal_nhit_layer_7 = ecal_nhit_layer_array[7];
	  b_ecal_nhit_layer_8 = ecal_nhit_layer_array[8];
	  b_ecal_nhit_layer_9 = ecal_nhit_layer_array[9];
	  b_ecal_nhit_layer_10 = ecal_nhit_layer_array[10];
	  b_ecal_nhit_layer_11 = ecal_nhit_layer_array[11];
	  b_ecal_nhit_layer_12 = ecal_nhit_layer_array[12];
	  b_ecal_nhit_layer_13 = ecal_nhit_layer_array[13];
	  b_ecal_nhit_layer_14 = ecal_nhit_layer_array[14];

	  b_ecal_nhit_layer_n_0 = ecal_nhit_layer_n_array[0];
	  b_ecal_nhit_layer_n_1 = ecal_nhit_layer_n_array[1];
	  b_ecal_nhit_layer_n_2 = ecal_nhit_layer_n_array[2];
	  b_ecal_nhit_layer_n_3 = ecal_nhit_layer_n_array[3];
	  b_ecal_nhit_layer_n_4 = ecal_nhit_layer_n_array[4];
	  b_ecal_nhit_layer_n_5 = ecal_nhit_layer_n_array[5];
	  b_ecal_nhit_layer_n_6 = ecal_nhit_layer_n_array[6];
	  b_ecal_nhit_layer_n_7 = ecal_nhit_layer_n_array[7];
	  b_ecal_nhit_layer_n_8 = ecal_nhit_layer_n_array[8];
	  b_ecal_nhit_layer_n_9 = ecal_nhit_layer_n_array[9];
	  b_ecal_nhit_layer_n_10 = ecal_nhit_layer_n_array[10];
	  b_ecal_nhit_layer_n_11 = ecal_nhit_layer_n_array[11];
	  b_ecal_nhit_layer_n_12 = ecal_nhit_layer_n_array[12];
	  b_ecal_nhit_layer_n_13 = ecal_nhit_layer_n_array[13];
	  b_ecal_nhit_layer_n_14 = ecal_nhit_layer_n_array[14];

	  b_ecal_shower_nhit_max_layer = ecal_nhit_ilayermax;
	  b_ecal_shower_nhit_start_layer = ecal_nhit_ilayerstart;
	  b_ecal_shower_nhit_end_layer = ecal_nhit_ilayerend;
	  b_ecal_shower_nhit_start_10_layer = ecal_nhit_ilayerstart_10;
	  b_ecal_shower_nhit_end_10_layer = ecal_nhit_ilayerend_10;
	  b_ecal_shower_nhit_average = ecal_nhit_shower_averagevalue;
	  b_ecal_shower_nhit_max = ecal_nhit_shower_maxvalue;

	  b_ecal_sume_layer_0 = ecal_sume_layer_array[0];
          b_ecal_sume_layer_1 = ecal_sume_layer_array[1];
          b_ecal_sume_layer_2 = ecal_sume_layer_array[2];
          b_ecal_sume_layer_3 = ecal_sume_layer_array[3];
          b_ecal_sume_layer_4 = ecal_sume_layer_array[4];
          b_ecal_sume_layer_5 = ecal_sume_layer_array[5];
          b_ecal_sume_layer_6 = ecal_sume_layer_array[6];
          b_ecal_sume_layer_7 = ecal_sume_layer_array[7];
          b_ecal_sume_layer_8 = ecal_sume_layer_array[8];
          b_ecal_sume_layer_9 = ecal_sume_layer_array[9];
          b_ecal_sume_layer_10 = ecal_sume_layer_array[10];
          b_ecal_sume_layer_11 = ecal_sume_layer_array[11];
          b_ecal_sume_layer_12 = ecal_sume_layer_array[12];
          b_ecal_sume_layer_13 = ecal_sume_layer_array[13];
          b_ecal_sume_layer_14 = ecal_sume_layer_array[14];

          b_ecal_sume_layer_n_0 = ecal_sume_layer_n_array[0];
          b_ecal_sume_layer_n_1 = ecal_sume_layer_n_array[1];
          b_ecal_sume_layer_n_2 = ecal_sume_layer_n_array[2];
          b_ecal_sume_layer_n_3 = ecal_sume_layer_n_array[3];
          b_ecal_sume_layer_n_4 = ecal_sume_layer_n_array[4];
          b_ecal_sume_layer_n_5 = ecal_sume_layer_n_array[5];
          b_ecal_sume_layer_n_6 = ecal_sume_layer_n_array[6];
          b_ecal_sume_layer_n_7 = ecal_sume_layer_n_array[7];
          b_ecal_sume_layer_n_8 = ecal_sume_layer_n_array[8];
          b_ecal_sume_layer_n_9 = ecal_sume_layer_n_array[9];
          b_ecal_sume_layer_n_10 = ecal_sume_layer_n_array[10];
          b_ecal_sume_layer_n_11 = ecal_sume_layer_n_array[11];
          b_ecal_sume_layer_n_12 = ecal_sume_layer_n_array[12];
          b_ecal_sume_layer_n_13 = ecal_sume_layer_n_array[13];
          b_ecal_sume_layer_n_14 = ecal_sume_layer_n_array[14];

          b_ecal_shower_sume_max_layer = ecal_sume_ilayermax;
          b_ecal_shower_sume_start_layer = ecal_sume_ilayerstart;
          b_ecal_shower_sume_end_layer = ecal_sume_ilayerend;
          b_ecal_shower_sume_start_10_layer = ecal_sume_ilayerstart_10;
          b_ecal_shower_sume_end_10_layer = ecal_sume_ilayerend_10;
          b_ecal_shower_sume_average = ecal_sume_shower_averagevalue;
          b_ecal_shower_sume_max = ecal_sume_shower_maxvalue;

	  b_ecal_weighte_layer_0 = ecal_weighte_layer_array[0];
	  b_ecal_weighte_layer_1 = ecal_weighte_layer_array[1];
	  b_ecal_weighte_layer_2 = ecal_weighte_layer_array[2];
	  b_ecal_weighte_layer_3 = ecal_weighte_layer_array[3];
	  b_ecal_weighte_layer_4 = ecal_weighte_layer_array[4];
	  b_ecal_weighte_layer_5 = ecal_weighte_layer_array[5];
	  b_ecal_weighte_layer_6 = ecal_weighte_layer_array[6];
	  b_ecal_weighte_layer_7 = ecal_weighte_layer_array[7];
	  b_ecal_weighte_layer_8 = ecal_weighte_layer_array[8];
	  b_ecal_weighte_layer_9 = ecal_weighte_layer_array[9];
	  b_ecal_weighte_layer_10 = ecal_weighte_layer_array[10];
	  b_ecal_weighte_layer_11 = ecal_weighte_layer_array[11];
	  b_ecal_weighte_layer_12 = ecal_weighte_layer_array[12];
	  b_ecal_weighte_layer_13 = ecal_weighte_layer_array[13];
	  b_ecal_weighte_layer_14 = ecal_weighte_layer_array[14];

	  b_ecal_weighte_layer_n_0 = ecal_weighte_layer_n_array[0];
	  b_ecal_weighte_layer_n_1 = ecal_weighte_layer_n_array[1];
	  b_ecal_weighte_layer_n_2 = ecal_weighte_layer_n_array[2];
	  b_ecal_weighte_layer_n_3 = ecal_weighte_layer_n_array[3];
	  b_ecal_weighte_layer_n_4 = ecal_weighte_layer_n_array[4];
	  b_ecal_weighte_layer_n_5 = ecal_weighte_layer_n_array[5];
	  b_ecal_weighte_layer_n_6 = ecal_weighte_layer_n_array[6];
	  b_ecal_weighte_layer_n_7 = ecal_weighte_layer_n_array[7];
	  b_ecal_weighte_layer_n_8 = ecal_weighte_layer_n_array[8];
	  b_ecal_weighte_layer_n_9 = ecal_weighte_layer_n_array[9];
	  b_ecal_weighte_layer_n_10 = ecal_weighte_layer_n_array[10];
	  b_ecal_weighte_layer_n_11 = ecal_weighte_layer_n_array[11];
	  b_ecal_weighte_layer_n_12 = ecal_weighte_layer_n_array[12];
	  b_ecal_weighte_layer_n_13 = ecal_weighte_layer_n_array[13];
	  b_ecal_weighte_layer_n_14 = ecal_weighte_layer_n_array[14];

	  b_ecal_shower_weighte_max_layer = ecal_weighte_ilayermax;
	  b_ecal_shower_weighte_start_layer = ecal_weighte_ilayerstart;
	  b_ecal_shower_weighte_end_layer = ecal_weighte_ilayerend;
	  b_ecal_shower_weighte_start_10_layer = ecal_weighte_ilayerstart_10;
	  b_ecal_shower_weighte_end_10_layer = ecal_weighte_ilayerend_10;
	  b_ecal_shower_weighte_average = ecal_weighte_shower_averagevalue;
	  b_ecal_shower_weighte_max = ecal_weighte_shower_maxvalue;

	  b_ecal_bar_x_layer_0 = ecal_bar_layer_array[0][0];
	  b_ecal_bar_x_layer_1 = ecal_bar_layer_array[1][0];
	  b_ecal_bar_x_layer_2 = ecal_bar_layer_array[2][0];
	  b_ecal_bar_x_layer_3 = ecal_bar_layer_array[3][0];
	  b_ecal_bar_x_layer_4 = ecal_bar_layer_array[4][0];
	  b_ecal_bar_x_layer_5 = ecal_bar_layer_array[5][0];
	  b_ecal_bar_x_layer_6 = ecal_bar_layer_array[6][0];
	  b_ecal_bar_x_layer_7 = ecal_bar_layer_array[7][0];
	  b_ecal_bar_x_layer_8 = ecal_bar_layer_array[8][0];
	  b_ecal_bar_x_layer_9 = ecal_bar_layer_array[9][0];
	  b_ecal_bar_x_layer_10 = ecal_bar_layer_array[10][0];
	  b_ecal_bar_x_layer_11 = ecal_bar_layer_array[11][0];
	  b_ecal_bar_x_layer_12 = ecal_bar_layer_array[12][0];
	  b_ecal_bar_x_layer_13 = ecal_bar_layer_array[13][0];
	  b_ecal_bar_x_layer_14 = ecal_bar_layer_array[14][0];

	  b_ecal_bar_y_layer_0 = ecal_bar_layer_array[0][1];
	  b_ecal_bar_y_layer_1 = ecal_bar_layer_array[1][1];
	  b_ecal_bar_y_layer_2 = ecal_bar_layer_array[2][1];
	  b_ecal_bar_y_layer_3 = ecal_bar_layer_array[3][1];
	  b_ecal_bar_y_layer_4 = ecal_bar_layer_array[4][1];
	  b_ecal_bar_y_layer_5 = ecal_bar_layer_array[5][1];
	  b_ecal_bar_y_layer_6 = ecal_bar_layer_array[6][1];
	  b_ecal_bar_y_layer_7 = ecal_bar_layer_array[7][1];
	  b_ecal_bar_y_layer_8 = ecal_bar_layer_array[8][1];
	  b_ecal_bar_y_layer_9 = ecal_bar_layer_array[9][1];
	  b_ecal_bar_y_layer_10 = ecal_bar_layer_array[10][1];
	  b_ecal_bar_y_layer_11 = ecal_bar_layer_array[11][1];
	  b_ecal_bar_y_layer_12 = ecal_bar_layer_array[12][1];
	  b_ecal_bar_y_layer_13 = ecal_bar_layer_array[13][1];
	  b_ecal_bar_y_layer_14 = ecal_bar_layer_array[14][1];

	  b_ecal_bar_r_layer_0 = ecal_bar_layer_array[0][2];
          b_ecal_bar_r_layer_1 = ecal_bar_layer_array[1][2];
          b_ecal_bar_r_layer_2 = ecal_bar_layer_array[2][2];
          b_ecal_bar_r_layer_3 = ecal_bar_layer_array[3][2];
          b_ecal_bar_r_layer_4 = ecal_bar_layer_array[4][2];
          b_ecal_bar_r_layer_5 = ecal_bar_layer_array[5][2];
          b_ecal_bar_r_layer_6 = ecal_bar_layer_array[6][2];
          b_ecal_bar_r_layer_7 = ecal_bar_layer_array[7][2];
          b_ecal_bar_r_layer_8 = ecal_bar_layer_array[8][2];
          b_ecal_bar_r_layer_9 = ecal_bar_layer_array[9][2];
          b_ecal_bar_r_layer_10 = ecal_bar_layer_array[10][2];
          b_ecal_bar_r_layer_11 = ecal_bar_layer_array[11][2];
          b_ecal_bar_r_layer_12 = ecal_bar_layer_array[12][2];
          b_ecal_bar_r_layer_13 = ecal_bar_layer_array[13][2];
          b_ecal_bar_r_layer_14 = ecal_bar_layer_array[14][2];

	  //cout<<"hcal nhit: "<<hcal_nhit<<endl;
	  b_hcal_nhit = hcal_nhit;
          b_hcal_sume = hcal_sume;
          b_hcal_weighte = hcal_weighte;
          b_hcal_bar_x = hcal_bar_xyzr[0];
          b_hcal_bar_y = hcal_bar_xyzr[1];
          b_hcal_bar_z = hcal_bar_xyzr[2];
          b_hcal_bar_r = hcal_bar_xyzr[3];
          b_hcal_mol = hcal_mol_value;
	  b_hcal_hits_max_distance = hcal_hits_max_distance_value;
          b_hcal_MIP_Likeness = hcal_MIP_Likeness_value;

          b_hcal_radius90_layer_0 = hcal_radius90_layer_array[0];
          b_hcal_radius90_layer_1 = hcal_radius90_layer_array[1];
          b_hcal_radius90_layer_2 = hcal_radius90_layer_array[2];
          b_hcal_radius90_layer_3 = hcal_radius90_layer_array[3];
          b_hcal_radius90_layer_4 = hcal_radius90_layer_array[4];
          b_hcal_radius90_layer_5 = hcal_radius90_layer_array[5];
          b_hcal_radius90_layer_6 = hcal_radius90_layer_array[6];
          b_hcal_radius90_layer_7 = hcal_radius90_layer_array[7];
          b_hcal_radius90_layer_8 = hcal_radius90_layer_array[8];
          b_hcal_radius90_layer_9 = hcal_radius90_layer_array[9];
          b_hcal_radius90_layer_10 = hcal_radius90_layer_array[10];
          b_hcal_radius90_layer_11 = hcal_radius90_layer_array[11];
          b_hcal_radius90_layer_12 = hcal_radius90_layer_array[12];
          b_hcal_radius90_layer_13 = hcal_radius90_layer_array[13];
          b_hcal_radius90_layer_14 = hcal_radius90_layer_array[14];

	  b_hcal_nhit_layer_0 = hcal_nhit_layer_array[0];
	  b_hcal_nhit_layer_1 = hcal_nhit_layer_array[1];
	  b_hcal_nhit_layer_2 = hcal_nhit_layer_array[2];
	  b_hcal_nhit_layer_3 = hcal_nhit_layer_array[3];
	  b_hcal_nhit_layer_4 = hcal_nhit_layer_array[4];
	  b_hcal_nhit_layer_5 = hcal_nhit_layer_array[5];
	  b_hcal_nhit_layer_6 = hcal_nhit_layer_array[6];
	  b_hcal_nhit_layer_7 = hcal_nhit_layer_array[7];
	  b_hcal_nhit_layer_8 = hcal_nhit_layer_array[8];
	  b_hcal_nhit_layer_9 = hcal_nhit_layer_array[9];
	  b_hcal_nhit_layer_10 = hcal_nhit_layer_array[10];
	  b_hcal_nhit_layer_11 = hcal_nhit_layer_array[11];
	  b_hcal_nhit_layer_12 = hcal_nhit_layer_array[12];
	  b_hcal_nhit_layer_13 = hcal_nhit_layer_array[13];
	  b_hcal_nhit_layer_14 = hcal_nhit_layer_array[14];
	  b_hcal_nhit_layer_15 = hcal_nhit_layer_array[15];
	  b_hcal_nhit_layer_16 = hcal_nhit_layer_array[16];
	  b_hcal_nhit_layer_17 = hcal_nhit_layer_array[17];
	  b_hcal_nhit_layer_18 = hcal_nhit_layer_array[18];
	  b_hcal_nhit_layer_19 = hcal_nhit_layer_array[19];
	  b_hcal_nhit_layer_20 = hcal_nhit_layer_array[20];
	  b_hcal_nhit_layer_21 = hcal_nhit_layer_array[21];
	  b_hcal_nhit_layer_22 = hcal_nhit_layer_array[22];
	  b_hcal_nhit_layer_23 = hcal_nhit_layer_array[23];
	  b_hcal_nhit_layer_24 = hcal_nhit_layer_array[24];
	  b_hcal_nhit_layer_25 = hcal_nhit_layer_array[25];
	  b_hcal_nhit_layer_26 = hcal_nhit_layer_array[26];
	  b_hcal_nhit_layer_27 = hcal_nhit_layer_array[27];
	  b_hcal_nhit_layer_28 = hcal_nhit_layer_array[28];
	  b_hcal_nhit_layer_29 = hcal_nhit_layer_array[29];
	  b_hcal_nhit_layer_30 = hcal_nhit_layer_array[30];
	  b_hcal_nhit_layer_31 = hcal_nhit_layer_array[31];
	  b_hcal_nhit_layer_32 = hcal_nhit_layer_array[32];
	  b_hcal_nhit_layer_33 = hcal_nhit_layer_array[33];
	  b_hcal_nhit_layer_34 = hcal_nhit_layer_array[34];
	  b_hcal_nhit_layer_35 = hcal_nhit_layer_array[35];
	  b_hcal_nhit_layer_36 = hcal_nhit_layer_array[36];
	  b_hcal_nhit_layer_37 = hcal_nhit_layer_array[37];
	  b_hcal_nhit_layer_38 = hcal_nhit_layer_array[38];
	  b_hcal_nhit_layer_39 = hcal_nhit_layer_array[39];
	  b_hcal_nhit_layer_40 = hcal_nhit_layer_array[40];
	  b_hcal_nhit_layer_n_0 = hcal_nhit_layer_n_array[0];
	  b_hcal_nhit_layer_n_1 = hcal_nhit_layer_n_array[1];
	  b_hcal_nhit_layer_n_2 = hcal_nhit_layer_n_array[2];
	  b_hcal_nhit_layer_n_3 = hcal_nhit_layer_n_array[3];
	  b_hcal_nhit_layer_n_4 = hcal_nhit_layer_n_array[4];
	  b_hcal_nhit_layer_n_5 = hcal_nhit_layer_n_array[5];
	  b_hcal_nhit_layer_n_6 = hcal_nhit_layer_n_array[6];
	  b_hcal_nhit_layer_n_7 = hcal_nhit_layer_n_array[7];
	  b_hcal_nhit_layer_n_8 = hcal_nhit_layer_n_array[8];
	  b_hcal_nhit_layer_n_9 = hcal_nhit_layer_n_array[9];
	  b_hcal_nhit_layer_n_10 = hcal_nhit_layer_n_array[10];
	  b_hcal_nhit_layer_n_11 = hcal_nhit_layer_n_array[11];
	  b_hcal_nhit_layer_n_12 = hcal_nhit_layer_n_array[12];
	  b_hcal_nhit_layer_n_13 = hcal_nhit_layer_n_array[13];
	  b_hcal_nhit_layer_n_14 = hcal_nhit_layer_n_array[14];
	  b_hcal_nhit_layer_n_15 = hcal_nhit_layer_n_array[15];
	  b_hcal_nhit_layer_n_16 = hcal_nhit_layer_n_array[16];
	  b_hcal_nhit_layer_n_17 = hcal_nhit_layer_n_array[17];
	  b_hcal_nhit_layer_n_18 = hcal_nhit_layer_n_array[18];
	  b_hcal_nhit_layer_n_19 = hcal_nhit_layer_n_array[19];
	  b_hcal_nhit_layer_n_20 = hcal_nhit_layer_n_array[20];
	  b_hcal_nhit_layer_n_21 = hcal_nhit_layer_n_array[21];
	  b_hcal_nhit_layer_n_22 = hcal_nhit_layer_n_array[22];
	  b_hcal_nhit_layer_n_23 = hcal_nhit_layer_n_array[23];
	  b_hcal_nhit_layer_n_24 = hcal_nhit_layer_n_array[24];
	  b_hcal_nhit_layer_n_25 = hcal_nhit_layer_n_array[25];
	  b_hcal_nhit_layer_n_26 = hcal_nhit_layer_n_array[26];
	  b_hcal_nhit_layer_n_27 = hcal_nhit_layer_n_array[27];
	  b_hcal_nhit_layer_n_28 = hcal_nhit_layer_n_array[28];
	  b_hcal_nhit_layer_n_29 = hcal_nhit_layer_n_array[29];
	  b_hcal_nhit_layer_n_30 = hcal_nhit_layer_n_array[30];
	  b_hcal_nhit_layer_n_31 = hcal_nhit_layer_n_array[31];
	  b_hcal_nhit_layer_n_32 = hcal_nhit_layer_n_array[32];
	  b_hcal_nhit_layer_n_33 = hcal_nhit_layer_n_array[33];
	  b_hcal_nhit_layer_n_34 = hcal_nhit_layer_n_array[34];
	  b_hcal_nhit_layer_n_35 = hcal_nhit_layer_n_array[35];
	  b_hcal_nhit_layer_n_36 = hcal_nhit_layer_n_array[36];
	  b_hcal_nhit_layer_n_37 = hcal_nhit_layer_n_array[37];
	  b_hcal_nhit_layer_n_38 = hcal_nhit_layer_n_array[38];
	  b_hcal_nhit_layer_n_39 = hcal_nhit_layer_n_array[39];
	  b_hcal_nhit_layer_n_40 = hcal_nhit_layer_n_array[40];

          b_hcal_shower_nhit_max_layer = hcal_nhit_ilayermax;
          b_hcal_shower_nhit_start_layer = hcal_nhit_ilayerstart;
          b_hcal_shower_nhit_end_layer = hcal_nhit_ilayerend;
          b_hcal_shower_nhit_start_10_layer = hcal_nhit_ilayerstart_10;
          b_hcal_shower_nhit_end_10_layer = hcal_nhit_ilayerend_10;
          b_hcal_shower_nhit_average = hcal_nhit_shower_averagevalue;
          b_hcal_shower_nhit_max = hcal_nhit_shower_maxvalue;

	  b_hcal_sume_layer_0 = hcal_sume_layer_array[0];
          b_hcal_sume_layer_1 = hcal_sume_layer_array[1];
          b_hcal_sume_layer_2 = hcal_sume_layer_array[2];
          b_hcal_sume_layer_3 = hcal_sume_layer_array[3];
          b_hcal_sume_layer_4 = hcal_sume_layer_array[4];
          b_hcal_sume_layer_5 = hcal_sume_layer_array[5];
          b_hcal_sume_layer_6 = hcal_sume_layer_array[6];
          b_hcal_sume_layer_7 = hcal_sume_layer_array[7];
          b_hcal_sume_layer_8 = hcal_sume_layer_array[8];
          b_hcal_sume_layer_9 = hcal_sume_layer_array[9];
          b_hcal_sume_layer_10 = hcal_sume_layer_array[10];
          b_hcal_sume_layer_11 = hcal_sume_layer_array[11];
          b_hcal_sume_layer_12 = hcal_sume_layer_array[12];
          b_hcal_sume_layer_13 = hcal_sume_layer_array[13];
          b_hcal_sume_layer_14 = hcal_sume_layer_array[14];
          b_hcal_sume_layer_15 = hcal_sume_layer_array[15];
          b_hcal_sume_layer_16 = hcal_sume_layer_array[16];
          b_hcal_sume_layer_17 = hcal_sume_layer_array[17];
          b_hcal_sume_layer_18 = hcal_sume_layer_array[18];
          b_hcal_sume_layer_19 = hcal_sume_layer_array[19];
          b_hcal_sume_layer_20 = hcal_sume_layer_array[20];
          b_hcal_sume_layer_21 = hcal_sume_layer_array[21];
          b_hcal_sume_layer_22 = hcal_sume_layer_array[22];
          b_hcal_sume_layer_23 = hcal_sume_layer_array[23];
          b_hcal_sume_layer_24 = hcal_sume_layer_array[24];
          b_hcal_sume_layer_25 = hcal_sume_layer_array[25];
          b_hcal_sume_layer_26 = hcal_sume_layer_array[26];
          b_hcal_sume_layer_27 = hcal_sume_layer_array[27];
          b_hcal_sume_layer_28 = hcal_sume_layer_array[28];
          b_hcal_sume_layer_29 = hcal_sume_layer_array[29];
          b_hcal_sume_layer_30 = hcal_sume_layer_array[30];
          b_hcal_sume_layer_31 = hcal_sume_layer_array[31];
          b_hcal_sume_layer_32 = hcal_sume_layer_array[32];
          b_hcal_sume_layer_33 = hcal_sume_layer_array[33];
          b_hcal_sume_layer_34 = hcal_sume_layer_array[34];
          b_hcal_sume_layer_35 = hcal_sume_layer_array[35];
          b_hcal_sume_layer_36 = hcal_sume_layer_array[36];
          b_hcal_sume_layer_37 = hcal_sume_layer_array[37];
          b_hcal_sume_layer_38 = hcal_sume_layer_array[38];
          b_hcal_sume_layer_39 = hcal_sume_layer_array[39];
          b_hcal_sume_layer_40 = hcal_sume_layer_array[40];
          b_hcal_sume_layer_n_0 = hcal_sume_layer_n_array[0];
          b_hcal_sume_layer_n_1 = hcal_sume_layer_n_array[1];
          b_hcal_sume_layer_n_2 = hcal_sume_layer_n_array[2];
          b_hcal_sume_layer_n_3 = hcal_sume_layer_n_array[3];
          b_hcal_sume_layer_n_4 = hcal_sume_layer_n_array[4];
          b_hcal_sume_layer_n_5 = hcal_sume_layer_n_array[5];
          b_hcal_sume_layer_n_6 = hcal_sume_layer_n_array[6];
          b_hcal_sume_layer_n_7 = hcal_sume_layer_n_array[7];
          b_hcal_sume_layer_n_8 = hcal_sume_layer_n_array[8];
          b_hcal_sume_layer_n_9 = hcal_sume_layer_n_array[9];
          b_hcal_sume_layer_n_10 = hcal_sume_layer_n_array[10];
          b_hcal_sume_layer_n_11 = hcal_sume_layer_n_array[11];
          b_hcal_sume_layer_n_12 = hcal_sume_layer_n_array[12];
          b_hcal_sume_layer_n_13 = hcal_sume_layer_n_array[13];
          b_hcal_sume_layer_n_14 = hcal_sume_layer_n_array[14];
          b_hcal_sume_layer_n_15 = hcal_sume_layer_n_array[15];
          b_hcal_sume_layer_n_16 = hcal_sume_layer_n_array[16];
          b_hcal_sume_layer_n_17 = hcal_sume_layer_n_array[17];
          b_hcal_sume_layer_n_18 = hcal_sume_layer_n_array[18];
          b_hcal_sume_layer_n_19 = hcal_sume_layer_n_array[19];
          b_hcal_sume_layer_n_20 = hcal_sume_layer_n_array[20];
          b_hcal_sume_layer_n_21 = hcal_sume_layer_n_array[21];
          b_hcal_sume_layer_n_22 = hcal_sume_layer_n_array[22];
          b_hcal_sume_layer_n_23 = hcal_sume_layer_n_array[23];
          b_hcal_sume_layer_n_24 = hcal_sume_layer_n_array[24];
          b_hcal_sume_layer_n_25 = hcal_sume_layer_n_array[25];
          b_hcal_sume_layer_n_26 = hcal_sume_layer_n_array[26];
          b_hcal_sume_layer_n_27 = hcal_sume_layer_n_array[27];
          b_hcal_sume_layer_n_28 = hcal_sume_layer_n_array[28];
          b_hcal_sume_layer_n_29 = hcal_sume_layer_n_array[29];
          b_hcal_sume_layer_n_30 = hcal_sume_layer_n_array[30];
          b_hcal_sume_layer_n_31 = hcal_sume_layer_n_array[31];
          b_hcal_sume_layer_n_32 = hcal_sume_layer_n_array[32];
          b_hcal_sume_layer_n_33 = hcal_sume_layer_n_array[33];
          b_hcal_sume_layer_n_34 = hcal_sume_layer_n_array[34];
          b_hcal_sume_layer_n_35 = hcal_sume_layer_n_array[35];
          b_hcal_sume_layer_n_36 = hcal_sume_layer_n_array[36];
          b_hcal_sume_layer_n_37 = hcal_sume_layer_n_array[37];
          b_hcal_sume_layer_n_38 = hcal_sume_layer_n_array[38];
          b_hcal_sume_layer_n_39 = hcal_sume_layer_n_array[39];
          b_hcal_sume_layer_n_40 = hcal_sume_layer_n_array[40];

          b_hcal_shower_sume_max_layer = hcal_sume_ilayermax;
          b_hcal_shower_sume_start_layer = hcal_sume_ilayerstart;
          b_hcal_shower_sume_end_layer = hcal_sume_ilayerend;
          b_hcal_shower_sume_start_10_layer = hcal_sume_ilayerstart_10;
          b_hcal_shower_sume_end_10_layer = hcal_sume_ilayerend_10;
          b_hcal_shower_sume_average = hcal_sume_shower_averagevalue;
          b_hcal_shower_sume_max = hcal_sume_shower_maxvalue;

	  b_hcal_weighte_layer_0 = hcal_weighte_layer_array[0];
          b_hcal_weighte_layer_1 = hcal_weighte_layer_array[1];
          b_hcal_weighte_layer_2 = hcal_weighte_layer_array[2];
          b_hcal_weighte_layer_3 = hcal_weighte_layer_array[3];
          b_hcal_weighte_layer_4 = hcal_weighte_layer_array[4];
          b_hcal_weighte_layer_5 = hcal_weighte_layer_array[5];
          b_hcal_weighte_layer_6 = hcal_weighte_layer_array[6];
          b_hcal_weighte_layer_7 = hcal_weighte_layer_array[7];
          b_hcal_weighte_layer_8 = hcal_weighte_layer_array[8];
          b_hcal_weighte_layer_9 = hcal_weighte_layer_array[9];
          b_hcal_weighte_layer_10 = hcal_weighte_layer_array[10];
          b_hcal_weighte_layer_11 = hcal_weighte_layer_array[11];
          b_hcal_weighte_layer_12 = hcal_weighte_layer_array[12];
          b_hcal_weighte_layer_13 = hcal_weighte_layer_array[13];
          b_hcal_weighte_layer_14 = hcal_weighte_layer_array[14];
          b_hcal_weighte_layer_15 = hcal_weighte_layer_array[15];
          b_hcal_weighte_layer_16 = hcal_weighte_layer_array[16];
          b_hcal_weighte_layer_17 = hcal_weighte_layer_array[17];
          b_hcal_weighte_layer_18 = hcal_weighte_layer_array[18];
          b_hcal_weighte_layer_19 = hcal_weighte_layer_array[19];
          b_hcal_weighte_layer_20 = hcal_weighte_layer_array[20];
          b_hcal_weighte_layer_21 = hcal_weighte_layer_array[21];
          b_hcal_weighte_layer_22 = hcal_weighte_layer_array[22];
          b_hcal_weighte_layer_23 = hcal_weighte_layer_array[23];
          b_hcal_weighte_layer_24 = hcal_weighte_layer_array[24];
          b_hcal_weighte_layer_25 = hcal_weighte_layer_array[25];
          b_hcal_weighte_layer_26 = hcal_weighte_layer_array[26];
          b_hcal_weighte_layer_27 = hcal_weighte_layer_array[27];
          b_hcal_weighte_layer_28 = hcal_weighte_layer_array[28];
          b_hcal_weighte_layer_29 = hcal_weighte_layer_array[29];
          b_hcal_weighte_layer_30 = hcal_weighte_layer_array[30];
          b_hcal_weighte_layer_31 = hcal_weighte_layer_array[31];
          b_hcal_weighte_layer_32 = hcal_weighte_layer_array[32];
          b_hcal_weighte_layer_33 = hcal_weighte_layer_array[33];
          b_hcal_weighte_layer_34 = hcal_weighte_layer_array[34];
          b_hcal_weighte_layer_35 = hcal_weighte_layer_array[35];
          b_hcal_weighte_layer_36 = hcal_weighte_layer_array[36];
          b_hcal_weighte_layer_37 = hcal_weighte_layer_array[37];
          b_hcal_weighte_layer_38 = hcal_weighte_layer_array[38];
          b_hcal_weighte_layer_39 = hcal_weighte_layer_array[39];
          b_hcal_weighte_layer_40 = hcal_weighte_layer_array[40];
          b_hcal_weighte_layer_n_0 = hcal_weighte_layer_n_array[0];
          b_hcal_weighte_layer_n_1 = hcal_weighte_layer_n_array[1];
          b_hcal_weighte_layer_n_2 = hcal_weighte_layer_n_array[2];
          b_hcal_weighte_layer_n_3 = hcal_weighte_layer_n_array[3];
          b_hcal_weighte_layer_n_4 = hcal_weighte_layer_n_array[4];
          b_hcal_weighte_layer_n_5 = hcal_weighte_layer_n_array[5];
          b_hcal_weighte_layer_n_6 = hcal_weighte_layer_n_array[6];
          b_hcal_weighte_layer_n_7 = hcal_weighte_layer_n_array[7];
          b_hcal_weighte_layer_n_8 = hcal_weighte_layer_n_array[8];
          b_hcal_weighte_layer_n_9 = hcal_weighte_layer_n_array[9];
          b_hcal_weighte_layer_n_10 = hcal_weighte_layer_n_array[10];
          b_hcal_weighte_layer_n_11 = hcal_weighte_layer_n_array[11];
          b_hcal_weighte_layer_n_12 = hcal_weighte_layer_n_array[12];
          b_hcal_weighte_layer_n_13 = hcal_weighte_layer_n_array[13];
          b_hcal_weighte_layer_n_14 = hcal_weighte_layer_n_array[14];
          b_hcal_weighte_layer_n_15 = hcal_weighte_layer_n_array[15];
          b_hcal_weighte_layer_n_16 = hcal_weighte_layer_n_array[16];
          b_hcal_weighte_layer_n_17 = hcal_weighte_layer_n_array[17];
          b_hcal_weighte_layer_n_18 = hcal_weighte_layer_n_array[18];
          b_hcal_weighte_layer_n_19 = hcal_weighte_layer_n_array[19];
          b_hcal_weighte_layer_n_20 = hcal_weighte_layer_n_array[20];
          b_hcal_weighte_layer_n_21 = hcal_weighte_layer_n_array[21];
          b_hcal_weighte_layer_n_22 = hcal_weighte_layer_n_array[22];
          b_hcal_weighte_layer_n_23 = hcal_weighte_layer_n_array[23];
          b_hcal_weighte_layer_n_24 = hcal_weighte_layer_n_array[24];
          b_hcal_weighte_layer_n_25 = hcal_weighte_layer_n_array[25];
          b_hcal_weighte_layer_n_26 = hcal_weighte_layer_n_array[26];
          b_hcal_weighte_layer_n_27 = hcal_weighte_layer_n_array[27];
          b_hcal_weighte_layer_n_28 = hcal_weighte_layer_n_array[28];
          b_hcal_weighte_layer_n_29 = hcal_weighte_layer_n_array[29];
          b_hcal_weighte_layer_n_30 = hcal_weighte_layer_n_array[30];
          b_hcal_weighte_layer_n_31 = hcal_weighte_layer_n_array[31];
          b_hcal_weighte_layer_n_32 = hcal_weighte_layer_n_array[32];
          b_hcal_weighte_layer_n_33 = hcal_weighte_layer_n_array[33];
          b_hcal_weighte_layer_n_34 = hcal_weighte_layer_n_array[34];
          b_hcal_weighte_layer_n_35 = hcal_weighte_layer_n_array[35];
          b_hcal_weighte_layer_n_36 = hcal_weighte_layer_n_array[36];
          b_hcal_weighte_layer_n_37 = hcal_weighte_layer_n_array[37];
          b_hcal_weighte_layer_n_38 = hcal_weighte_layer_n_array[38];
          b_hcal_weighte_layer_n_39 = hcal_weighte_layer_n_array[39];
          b_hcal_weighte_layer_n_40 = hcal_weighte_layer_n_array[40];

          b_hcal_shower_weighte_max_layer = hcal_weighte_ilayermax;
          b_hcal_shower_weighte_start_layer = hcal_weighte_ilayerstart;
          b_hcal_shower_weighte_end_layer = hcal_weighte_ilayerend;
          b_hcal_shower_weighte_start_10_layer = hcal_weighte_ilayerstart_10;
          b_hcal_shower_weighte_end_10_layer = hcal_weighte_ilayerend_10;
          b_hcal_shower_weighte_average = hcal_weighte_shower_averagevalue;
          b_hcal_shower_weighte_max = hcal_weighte_shower_maxvalue;

          b_hcal_bar_x_layer_0 = hcal_bar_layer_array[0][0];
          b_hcal_bar_x_layer_1 = hcal_bar_layer_array[1][0];
          b_hcal_bar_x_layer_2 = hcal_bar_layer_array[2][0];
          b_hcal_bar_x_layer_3 = hcal_bar_layer_array[3][0];
          b_hcal_bar_x_layer_4 = hcal_bar_layer_array[4][0];
          b_hcal_bar_x_layer_5 = hcal_bar_layer_array[5][0];
          b_hcal_bar_x_layer_6 = hcal_bar_layer_array[6][0];
          b_hcal_bar_x_layer_7 = hcal_bar_layer_array[7][0];
          b_hcal_bar_x_layer_8 = hcal_bar_layer_array[8][0];
          b_hcal_bar_x_layer_9 = hcal_bar_layer_array[9][0];
          b_hcal_bar_x_layer_10 = hcal_bar_layer_array[10][0];
          b_hcal_bar_x_layer_11 = hcal_bar_layer_array[11][0];
          b_hcal_bar_x_layer_12 = hcal_bar_layer_array[12][0];
          b_hcal_bar_x_layer_13 = hcal_bar_layer_array[13][0];
          b_hcal_bar_x_layer_14 = hcal_bar_layer_array[14][0];

          b_hcal_bar_y_layer_0 = hcal_bar_layer_array[0][1];
          b_hcal_bar_y_layer_1 = hcal_bar_layer_array[1][1];
          b_hcal_bar_y_layer_2 = hcal_bar_layer_array[2][1];
          b_hcal_bar_y_layer_3 = hcal_bar_layer_array[3][1];
          b_hcal_bar_y_layer_4 = hcal_bar_layer_array[4][1];
          b_hcal_bar_y_layer_5 = hcal_bar_layer_array[5][1];
          b_hcal_bar_y_layer_6 = hcal_bar_layer_array[6][1];
          b_hcal_bar_y_layer_7 = hcal_bar_layer_array[7][1];
          b_hcal_bar_y_layer_8 = hcal_bar_layer_array[8][1];
          b_hcal_bar_y_layer_9 = hcal_bar_layer_array[9][1];
          b_hcal_bar_y_layer_10 = hcal_bar_layer_array[10][1];
          b_hcal_bar_y_layer_11 = hcal_bar_layer_array[11][1];
          b_hcal_bar_y_layer_12 = hcal_bar_layer_array[12][1];
          b_hcal_bar_y_layer_13 = hcal_bar_layer_array[13][1];
          b_hcal_bar_y_layer_14 = hcal_bar_layer_array[14][1];

          b_hcal_bar_r_layer_0 = hcal_bar_layer_array[0][2];
          b_hcal_bar_r_layer_1 = hcal_bar_layer_array[1][2];
          b_hcal_bar_r_layer_2 = hcal_bar_layer_array[2][2];
          b_hcal_bar_r_layer_3 = hcal_bar_layer_array[3][2];
          b_hcal_bar_r_layer_4 = hcal_bar_layer_array[4][2];
          b_hcal_bar_r_layer_5 = hcal_bar_layer_array[5][2];
          b_hcal_bar_r_layer_6 = hcal_bar_layer_array[6][2];
          b_hcal_bar_r_layer_7 = hcal_bar_layer_array[7][2];
          b_hcal_bar_r_layer_8 = hcal_bar_layer_array[8][2];
          b_hcal_bar_r_layer_9 = hcal_bar_layer_array[9][2];
          b_hcal_bar_r_layer_10 = hcal_bar_layer_array[10][2];
          b_hcal_bar_r_layer_11 = hcal_bar_layer_array[11][2];
          b_hcal_bar_r_layer_12 = hcal_bar_layer_array[12][2];
	  
	  
	  b_total_nhit = total_nhit;
          b_total_sume = total_sume;
          b_total_weighte = total_weighte;
          b_total_bar_x = total_bar_xyzr[0];
          b_total_bar_y = total_bar_xyzr[1];
          b_total_bar_z = total_bar_xyzr[2];
          b_total_bar_r = total_bar_xyzr[3];
          b_total_mol = total_mol_value;
          b_total_hits_max_distance = total_hits_max_distance_value;
	  b_total_MIP_Likeness = total_MIP_Likeness_value;

	  b_total_radius90_layer_0 = total_radius90_layer_array[0];
          b_total_radius90_layer_1 = total_radius90_layer_array[1];
          b_total_radius90_layer_2 = total_radius90_layer_array[2];
          b_total_radius90_layer_3 = total_radius90_layer_array[3];
          b_total_radius90_layer_4 = total_radius90_layer_array[4];
          b_total_radius90_layer_5 = total_radius90_layer_array[5];
          b_total_radius90_layer_6 = total_radius90_layer_array[6];
          b_total_radius90_layer_7 = total_radius90_layer_array[7];
          b_total_radius90_layer_8 = total_radius90_layer_array[8];
          b_total_radius90_layer_9 = total_radius90_layer_array[9];
          b_total_radius90_layer_10 = total_radius90_layer_array[10];
          b_total_radius90_layer_11 = total_radius90_layer_array[11];
          b_total_radius90_layer_12 = total_radius90_layer_array[12];
          b_total_radius90_layer_13 = total_radius90_layer_array[13];
          b_total_radius90_layer_14 = total_radius90_layer_array[14];

	  b_total_nhit_layer_0 = total_nhit_layer_array[0];
	  b_total_nhit_layer_1 = total_nhit_layer_array[1];
	  b_total_nhit_layer_2 = total_nhit_layer_array[2];
	  b_total_nhit_layer_3 = total_nhit_layer_array[3];
	  b_total_nhit_layer_4 = total_nhit_layer_array[4];
	  b_total_nhit_layer_5 = total_nhit_layer_array[5];
	  b_total_nhit_layer_6 = total_nhit_layer_array[6];
	  b_total_nhit_layer_7 = total_nhit_layer_array[7];
	  b_total_nhit_layer_8 = total_nhit_layer_array[8];
	  b_total_nhit_layer_9 = total_nhit_layer_array[9];
	  b_total_nhit_layer_10 = total_nhit_layer_array[10];
	  b_total_nhit_layer_11 = total_nhit_layer_array[11];
	  b_total_nhit_layer_12 = total_nhit_layer_array[12];
	  b_total_nhit_layer_13 = total_nhit_layer_array[13];
	  b_total_nhit_layer_14 = total_nhit_layer_array[14];
	  b_total_nhit_layer_15 = total_nhit_layer_array[15];
	  b_total_nhit_layer_16 = total_nhit_layer_array[16];
	  b_total_nhit_layer_17 = total_nhit_layer_array[17];
	  b_total_nhit_layer_18 = total_nhit_layer_array[18];
	  b_total_nhit_layer_19 = total_nhit_layer_array[19];
	  b_total_nhit_layer_20 = total_nhit_layer_array[20];
	  b_total_nhit_layer_21 = total_nhit_layer_array[21];
	  b_total_nhit_layer_22 = total_nhit_layer_array[22];
	  b_total_nhit_layer_23 = total_nhit_layer_array[23];
	  b_total_nhit_layer_24 = total_nhit_layer_array[24];
	  b_total_nhit_layer_25 = total_nhit_layer_array[25];
	  b_total_nhit_layer_26 = total_nhit_layer_array[26];
	  b_total_nhit_layer_27 = total_nhit_layer_array[27];
	  b_total_nhit_layer_28 = total_nhit_layer_array[28];
	  b_total_nhit_layer_29 = total_nhit_layer_array[29];
	  b_total_nhit_layer_30 = total_nhit_layer_array[30];
	  b_total_nhit_layer_31 = total_nhit_layer_array[31];
	  b_total_nhit_layer_32 = total_nhit_layer_array[32];
	  b_total_nhit_layer_33 = total_nhit_layer_array[33];
	  b_total_nhit_layer_34 = total_nhit_layer_array[34];
	  b_total_nhit_layer_35 = total_nhit_layer_array[35];
	  b_total_nhit_layer_36 = total_nhit_layer_array[36];
	  b_total_nhit_layer_37 = total_nhit_layer_array[37];
	  b_total_nhit_layer_38 = total_nhit_layer_array[38];
	  b_total_nhit_layer_39 = total_nhit_layer_array[39];
	  b_total_nhit_layer_40 = total_nhit_layer_array[40];
	  b_total_nhit_layer_41 = total_nhit_layer_array[41];
	  b_total_nhit_layer_42 = total_nhit_layer_array[42];
	  b_total_nhit_layer_43 = total_nhit_layer_array[43];
	  b_total_nhit_layer_44 = total_nhit_layer_array[44];
	  b_total_nhit_layer_45 = total_nhit_layer_array[45];
	  b_total_nhit_layer_46 = total_nhit_layer_array[46];
	  b_total_nhit_layer_47 = total_nhit_layer_array[47];
	  b_total_nhit_layer_48 = total_nhit_layer_array[48];
	  b_total_nhit_layer_49 = total_nhit_layer_array[49];
	  b_total_nhit_layer_50 = total_nhit_layer_array[50];
	  b_total_nhit_layer_51 = total_nhit_layer_array[51];
	  b_total_nhit_layer_52 = total_nhit_layer_array[52];
	  b_total_nhit_layer_53 = total_nhit_layer_array[53];
	  b_total_nhit_layer_54 = total_nhit_layer_array[54];
	  b_total_nhit_layer_55 = total_nhit_layer_array[55];
	  b_total_nhit_layer_n_0 = total_nhit_layer_n_array[0];
	  b_total_nhit_layer_n_1 = total_nhit_layer_n_array[1];
	  b_total_nhit_layer_n_2 = total_nhit_layer_n_array[2];
	  b_total_nhit_layer_n_3 = total_nhit_layer_n_array[3];
	  b_total_nhit_layer_n_4 = total_nhit_layer_n_array[4];
	  b_total_nhit_layer_n_5 = total_nhit_layer_n_array[5];
	  b_total_nhit_layer_n_6 = total_nhit_layer_n_array[6];
	  b_total_nhit_layer_n_7 = total_nhit_layer_n_array[7];
	  b_total_nhit_layer_n_8 = total_nhit_layer_n_array[8];
	  b_total_nhit_layer_n_9 = total_nhit_layer_n_array[9];
	  b_total_nhit_layer_n_10 = total_nhit_layer_n_array[10];
	  b_total_nhit_layer_n_11 = total_nhit_layer_n_array[11];
	  b_total_nhit_layer_n_12 = total_nhit_layer_n_array[12];
	  b_total_nhit_layer_n_13 = total_nhit_layer_n_array[13];
	  b_total_nhit_layer_n_14 = total_nhit_layer_n_array[14];
	  b_total_nhit_layer_n_15 = total_nhit_layer_n_array[15];
	  b_total_nhit_layer_n_16 = total_nhit_layer_n_array[16];
	  b_total_nhit_layer_n_17 = total_nhit_layer_n_array[17];
	  b_total_nhit_layer_n_18 = total_nhit_layer_n_array[18];
	  b_total_nhit_layer_n_19 = total_nhit_layer_n_array[19];
	  b_total_nhit_layer_n_20 = total_nhit_layer_n_array[20];
	  b_total_nhit_layer_n_21 = total_nhit_layer_n_array[21];
	  b_total_nhit_layer_n_22 = total_nhit_layer_n_array[22];
	  b_total_nhit_layer_n_23 = total_nhit_layer_n_array[23];
	  b_total_nhit_layer_n_24 = total_nhit_layer_n_array[24];
	  b_total_nhit_layer_n_25 = total_nhit_layer_n_array[25];
	  b_total_nhit_layer_n_26 = total_nhit_layer_n_array[26];
	  b_total_nhit_layer_n_27 = total_nhit_layer_n_array[27];
	  b_total_nhit_layer_n_28 = total_nhit_layer_n_array[28];
	  b_total_nhit_layer_n_29 = total_nhit_layer_n_array[29];
	  b_total_nhit_layer_n_30 = total_nhit_layer_n_array[30];
	  b_total_nhit_layer_n_31 = total_nhit_layer_n_array[31];
	  b_total_nhit_layer_n_32 = total_nhit_layer_n_array[32];
	  b_total_nhit_layer_n_33 = total_nhit_layer_n_array[33];
	  b_total_nhit_layer_n_34 = total_nhit_layer_n_array[34];
	  b_total_nhit_layer_n_35 = total_nhit_layer_n_array[35];
	  b_total_nhit_layer_n_36 = total_nhit_layer_n_array[36];
	  b_total_nhit_layer_n_37 = total_nhit_layer_n_array[37];
	  b_total_nhit_layer_n_38 = total_nhit_layer_n_array[38];
	  b_total_nhit_layer_n_39 = total_nhit_layer_n_array[39];
	  b_total_nhit_layer_n_40 = total_nhit_layer_n_array[40];
	  b_total_nhit_layer_n_41 = total_nhit_layer_n_array[41];
	  b_total_nhit_layer_n_42 = total_nhit_layer_n_array[42];
	  b_total_nhit_layer_n_43 = total_nhit_layer_n_array[43];
	  b_total_nhit_layer_n_44 = total_nhit_layer_n_array[44];
	  b_total_nhit_layer_n_45 = total_nhit_layer_n_array[45];
	  b_total_nhit_layer_n_46 = total_nhit_layer_n_array[46];
	  b_total_nhit_layer_n_47 = total_nhit_layer_n_array[47];
	  b_total_nhit_layer_n_48 = total_nhit_layer_n_array[48];
	  b_total_nhit_layer_n_49 = total_nhit_layer_n_array[49];
	  b_total_nhit_layer_n_50 = total_nhit_layer_n_array[50];
	  b_total_nhit_layer_n_51 = total_nhit_layer_n_array[51];
	  b_total_nhit_layer_n_52 = total_nhit_layer_n_array[52];
	  b_total_nhit_layer_n_53 = total_nhit_layer_n_array[53];
	  b_total_nhit_layer_n_54 = total_nhit_layer_n_array[54];
	  b_total_nhit_layer_n_55 = total_nhit_layer_n_array[55];

          b_total_shower_nhit_max_layer = total_nhit_ilayermax;
          b_total_shower_nhit_start_layer = total_nhit_ilayerstart;
          b_total_shower_nhit_end_layer = total_nhit_ilayerend;
          b_total_shower_nhit_start_10_layer = total_nhit_ilayerstart_10;
          b_total_shower_nhit_end_10_layer = total_nhit_ilayerend_10;
          b_total_shower_nhit_average = total_nhit_shower_averagevalue;
          b_total_shower_nhit_max = total_nhit_shower_maxvalue;

	  b_total_sume_layer_0 = total_sume_layer_array[0];
          b_total_sume_layer_1 = total_sume_layer_array[1];
          b_total_sume_layer_2 = total_sume_layer_array[2];
          b_total_sume_layer_3 = total_sume_layer_array[3];
          b_total_sume_layer_4 = total_sume_layer_array[4];
          b_total_sume_layer_5 = total_sume_layer_array[5];
          b_total_sume_layer_6 = total_sume_layer_array[6];
          b_total_sume_layer_7 = total_sume_layer_array[7];
          b_total_sume_layer_8 = total_sume_layer_array[8];
          b_total_sume_layer_9 = total_sume_layer_array[9];
          b_total_sume_layer_10 = total_sume_layer_array[10];
          b_total_sume_layer_11 = total_sume_layer_array[11];
          b_total_sume_layer_12 = total_sume_layer_array[12];
          b_total_sume_layer_13 = total_sume_layer_array[13];
          b_total_sume_layer_14 = total_sume_layer_array[14];
          b_total_sume_layer_15 = total_sume_layer_array[15];
          b_total_sume_layer_16 = total_sume_layer_array[16];
          b_total_sume_layer_17 = total_sume_layer_array[17];
          b_total_sume_layer_18 = total_sume_layer_array[18];
          b_total_sume_layer_19 = total_sume_layer_array[19];
          b_total_sume_layer_20 = total_sume_layer_array[20];
          b_total_sume_layer_21 = total_sume_layer_array[21];
          b_total_sume_layer_22 = total_sume_layer_array[22];
          b_total_sume_layer_23 = total_sume_layer_array[23];
          b_total_sume_layer_24 = total_sume_layer_array[24];
          b_total_sume_layer_25 = total_sume_layer_array[25];
          b_total_sume_layer_26 = total_sume_layer_array[26];
          b_total_sume_layer_27 = total_sume_layer_array[27];
          b_total_sume_layer_28 = total_sume_layer_array[28];
          b_total_sume_layer_29 = total_sume_layer_array[29];
          b_total_sume_layer_30 = total_sume_layer_array[30];
          b_total_sume_layer_31 = total_sume_layer_array[31];
          b_total_sume_layer_32 = total_sume_layer_array[32];
          b_total_sume_layer_33 = total_sume_layer_array[33];
          b_total_sume_layer_34 = total_sume_layer_array[34];
          b_total_sume_layer_35 = total_sume_layer_array[35];
          b_total_sume_layer_36 = total_sume_layer_array[36];
          b_total_sume_layer_37 = total_sume_layer_array[37];
          b_total_sume_layer_38 = total_sume_layer_array[38];
          b_total_sume_layer_39 = total_sume_layer_array[39];
          b_total_sume_layer_40 = total_sume_layer_array[40];
          b_total_sume_layer_41 = total_sume_layer_array[41];
          b_total_sume_layer_42 = total_sume_layer_array[42];
          b_total_sume_layer_43 = total_sume_layer_array[43];
          b_total_sume_layer_44 = total_sume_layer_array[44];
          b_total_sume_layer_45 = total_sume_layer_array[45];
          b_total_sume_layer_46 = total_sume_layer_array[46];
          b_total_sume_layer_47 = total_sume_layer_array[47];
          b_total_sume_layer_48 = total_sume_layer_array[48];
          b_total_sume_layer_49 = total_sume_layer_array[49];
          b_total_sume_layer_50 = total_sume_layer_array[50];
          b_total_sume_layer_51 = total_sume_layer_array[51];
          b_total_sume_layer_52 = total_sume_layer_array[52];
          b_total_sume_layer_53 = total_sume_layer_array[53];
          b_total_sume_layer_54 = total_sume_layer_array[54];
          b_total_sume_layer_55 = total_sume_layer_array[55];
          b_total_sume_layer_n_0 = total_sume_layer_n_array[0];
          b_total_sume_layer_n_1 = total_sume_layer_n_array[1];
          b_total_sume_layer_n_2 = total_sume_layer_n_array[2];
          b_total_sume_layer_n_3 = total_sume_layer_n_array[3];
          b_total_sume_layer_n_4 = total_sume_layer_n_array[4];
          b_total_sume_layer_n_5 = total_sume_layer_n_array[5];
          b_total_sume_layer_n_6 = total_sume_layer_n_array[6];
          b_total_sume_layer_n_7 = total_sume_layer_n_array[7];
          b_total_sume_layer_n_8 = total_sume_layer_n_array[8];
          b_total_sume_layer_n_9 = total_sume_layer_n_array[9];
          b_total_sume_layer_n_10 = total_sume_layer_n_array[10];
          b_total_sume_layer_n_11 = total_sume_layer_n_array[11];
          b_total_sume_layer_n_12 = total_sume_layer_n_array[12];
          b_total_sume_layer_n_13 = total_sume_layer_n_array[13];
          b_total_sume_layer_n_14 = total_sume_layer_n_array[14];
          b_total_sume_layer_n_15 = total_sume_layer_n_array[15];
          b_total_sume_layer_n_16 = total_sume_layer_n_array[16];
          b_total_sume_layer_n_17 = total_sume_layer_n_array[17];
          b_total_sume_layer_n_18 = total_sume_layer_n_array[18];
          b_total_sume_layer_n_19 = total_sume_layer_n_array[19];
          b_total_sume_layer_n_20 = total_sume_layer_n_array[20];
          b_total_sume_layer_n_21 = total_sume_layer_n_array[21];
          b_total_sume_layer_n_22 = total_sume_layer_n_array[22];
          b_total_sume_layer_n_23 = total_sume_layer_n_array[23];
          b_total_sume_layer_n_24 = total_sume_layer_n_array[24];
          b_total_sume_layer_n_25 = total_sume_layer_n_array[25];
          b_total_sume_layer_n_26 = total_sume_layer_n_array[26];
          b_total_sume_layer_n_27 = total_sume_layer_n_array[27];
          b_total_sume_layer_n_28 = total_sume_layer_n_array[28];
          b_total_sume_layer_n_29 = total_sume_layer_n_array[29];
          b_total_sume_layer_n_30 = total_sume_layer_n_array[30];
          b_total_sume_layer_n_31 = total_sume_layer_n_array[31];
          b_total_sume_layer_n_32 = total_sume_layer_n_array[32];
          b_total_sume_layer_n_33 = total_sume_layer_n_array[33];
          b_total_sume_layer_n_34 = total_sume_layer_n_array[34];
          b_total_sume_layer_n_35 = total_sume_layer_n_array[35];
          b_total_sume_layer_n_36 = total_sume_layer_n_array[36];
          b_total_sume_layer_n_37 = total_sume_layer_n_array[37];
          b_total_sume_layer_n_38 = total_sume_layer_n_array[38];
          b_total_sume_layer_n_39 = total_sume_layer_n_array[39];
          b_total_sume_layer_n_40 = total_sume_layer_n_array[40];
          b_total_sume_layer_n_41 = total_sume_layer_n_array[41];
          b_total_sume_layer_n_42 = total_sume_layer_n_array[42];
          b_total_sume_layer_n_43 = total_sume_layer_n_array[43];
          b_total_sume_layer_n_44 = total_sume_layer_n_array[44];
          b_total_sume_layer_n_45 = total_sume_layer_n_array[45];
          b_total_sume_layer_n_46 = total_sume_layer_n_array[46];
          b_total_sume_layer_n_47 = total_sume_layer_n_array[47];
          b_total_sume_layer_n_48 = total_sume_layer_n_array[48];
          b_total_sume_layer_n_49 = total_sume_layer_n_array[49];
          b_total_sume_layer_n_50 = total_sume_layer_n_array[50];
          b_total_sume_layer_n_51 = total_sume_layer_n_array[51];
          b_total_sume_layer_n_52 = total_sume_layer_n_array[52];
          b_total_sume_layer_n_53 = total_sume_layer_n_array[53];
          b_total_sume_layer_n_54 = total_sume_layer_n_array[54];
          b_total_sume_layer_n_55 = total_sume_layer_n_array[55];
	  
          b_total_shower_sume_max_layer = total_sume_ilayermax;
          b_total_shower_sume_start_layer = total_sume_ilayerstart;
          b_total_shower_sume_end_layer = total_sume_ilayerend;
          b_total_shower_sume_start_10_layer = total_sume_ilayerstart_10;
          b_total_shower_sume_end_10_layer = total_sume_ilayerend_10;
          b_total_shower_sume_average = total_sume_shower_averagevalue;
          b_total_shower_sume_max = total_sume_shower_maxvalue;

	  b_total_weighte_layer_0 = total_weighte_layer_array[0];
          b_total_weighte_layer_1 = total_weighte_layer_array[1];
          b_total_weighte_layer_2 = total_weighte_layer_array[2];
          b_total_weighte_layer_3 = total_weighte_layer_array[3];
          b_total_weighte_layer_4 = total_weighte_layer_array[4];
          b_total_weighte_layer_5 = total_weighte_layer_array[5];
          b_total_weighte_layer_6 = total_weighte_layer_array[6];
          b_total_weighte_layer_7 = total_weighte_layer_array[7];
          b_total_weighte_layer_8 = total_weighte_layer_array[8];
          b_total_weighte_layer_9 = total_weighte_layer_array[9];
          b_total_weighte_layer_10 = total_weighte_layer_array[10];
          b_total_weighte_layer_11 = total_weighte_layer_array[11];
          b_total_weighte_layer_12 = total_weighte_layer_array[12];
          b_total_weighte_layer_13 = total_weighte_layer_array[13];
          b_total_weighte_layer_14 = total_weighte_layer_array[14];
          b_total_weighte_layer_15 = total_weighte_layer_array[15];
          b_total_weighte_layer_16 = total_weighte_layer_array[16];
          b_total_weighte_layer_17 = total_weighte_layer_array[17];
          b_total_weighte_layer_18 = total_weighte_layer_array[18];
          b_total_weighte_layer_19 = total_weighte_layer_array[19];
          b_total_weighte_layer_20 = total_weighte_layer_array[20];
          b_total_weighte_layer_21 = total_weighte_layer_array[21];
          b_total_weighte_layer_22 = total_weighte_layer_array[22];
          b_total_weighte_layer_23 = total_weighte_layer_array[23];
          b_total_weighte_layer_24 = total_weighte_layer_array[24];
          b_total_weighte_layer_25 = total_weighte_layer_array[25];
          b_total_weighte_layer_26 = total_weighte_layer_array[26];
          b_total_weighte_layer_27 = total_weighte_layer_array[27];
          b_total_weighte_layer_28 = total_weighte_layer_array[28];
          b_total_weighte_layer_29 = total_weighte_layer_array[29];
          b_total_weighte_layer_30 = total_weighte_layer_array[30];
          b_total_weighte_layer_31 = total_weighte_layer_array[31];
          b_total_weighte_layer_32 = total_weighte_layer_array[32];
          b_total_weighte_layer_33 = total_weighte_layer_array[33];
          b_total_weighte_layer_34 = total_weighte_layer_array[34];
          b_total_weighte_layer_35 = total_weighte_layer_array[35];
          b_total_weighte_layer_36 = total_weighte_layer_array[36];
          b_total_weighte_layer_37 = total_weighte_layer_array[37];
          b_total_weighte_layer_38 = total_weighte_layer_array[38];
          b_total_weighte_layer_39 = total_weighte_layer_array[39];
          b_total_weighte_layer_40 = total_weighte_layer_array[40];
          b_total_weighte_layer_41 = total_weighte_layer_array[41];
          b_total_weighte_layer_42 = total_weighte_layer_array[42];
          b_total_weighte_layer_43 = total_weighte_layer_array[43];
          b_total_weighte_layer_44 = total_weighte_layer_array[44];
          b_total_weighte_layer_45 = total_weighte_layer_array[45];
          b_total_weighte_layer_46 = total_weighte_layer_array[46];
          b_total_weighte_layer_47 = total_weighte_layer_array[47];
          b_total_weighte_layer_48 = total_weighte_layer_array[48];
          b_total_weighte_layer_49 = total_weighte_layer_array[49];
          b_total_weighte_layer_50 = total_weighte_layer_array[50];
          b_total_weighte_layer_51 = total_weighte_layer_array[51];
          b_total_weighte_layer_52 = total_weighte_layer_array[52];
          b_total_weighte_layer_53 = total_weighte_layer_array[53];
          b_total_weighte_layer_54 = total_weighte_layer_array[54];
          b_total_weighte_layer_55 = total_weighte_layer_array[55];
          b_total_weighte_layer_n_0 = total_weighte_layer_n_array[0];
          b_total_weighte_layer_n_1 = total_weighte_layer_n_array[1];
          b_total_weighte_layer_n_2 = total_weighte_layer_n_array[2];
          b_total_weighte_layer_n_3 = total_weighte_layer_n_array[3];
          b_total_weighte_layer_n_4 = total_weighte_layer_n_array[4];
          b_total_weighte_layer_n_5 = total_weighte_layer_n_array[5];
          b_total_weighte_layer_n_6 = total_weighte_layer_n_array[6];
          b_total_weighte_layer_n_7 = total_weighte_layer_n_array[7];
          b_total_weighte_layer_n_8 = total_weighte_layer_n_array[8];
          b_total_weighte_layer_n_9 = total_weighte_layer_n_array[9];
          b_total_weighte_layer_n_10 = total_weighte_layer_n_array[10];
          b_total_weighte_layer_n_11 = total_weighte_layer_n_array[11];
          b_total_weighte_layer_n_12 = total_weighte_layer_n_array[12];
          b_total_weighte_layer_n_13 = total_weighte_layer_n_array[13];
          b_total_weighte_layer_n_14 = total_weighte_layer_n_array[14];
          b_total_weighte_layer_n_15 = total_weighte_layer_n_array[15];
          b_total_weighte_layer_n_16 = total_weighte_layer_n_array[16];
          b_total_weighte_layer_n_17 = total_weighte_layer_n_array[17];
          b_total_weighte_layer_n_18 = total_weighte_layer_n_array[18];
          b_total_weighte_layer_n_19 = total_weighte_layer_n_array[19];
          b_total_weighte_layer_n_20 = total_weighte_layer_n_array[20];
          b_total_weighte_layer_n_21 = total_weighte_layer_n_array[21];
          b_total_weighte_layer_n_22 = total_weighte_layer_n_array[22];
          b_total_weighte_layer_n_23 = total_weighte_layer_n_array[23];
          b_total_weighte_layer_n_24 = total_weighte_layer_n_array[24];
          b_total_weighte_layer_n_25 = total_weighte_layer_n_array[25];
          b_total_weighte_layer_n_26 = total_weighte_layer_n_array[26];
          b_total_weighte_layer_n_27 = total_weighte_layer_n_array[27];
          b_total_weighte_layer_n_28 = total_weighte_layer_n_array[28];
          b_total_weighte_layer_n_29 = total_weighte_layer_n_array[29];
          b_total_weighte_layer_n_30 = total_weighte_layer_n_array[30];
          b_total_weighte_layer_n_31 = total_weighte_layer_n_array[31];
          b_total_weighte_layer_n_32 = total_weighte_layer_n_array[32];
          b_total_weighte_layer_n_33 = total_weighte_layer_n_array[33];
          b_total_weighte_layer_n_34 = total_weighte_layer_n_array[34];
          b_total_weighte_layer_n_35 = total_weighte_layer_n_array[35];
          b_total_weighte_layer_n_36 = total_weighte_layer_n_array[36];
          b_total_weighte_layer_n_37 = total_weighte_layer_n_array[37];
          b_total_weighte_layer_n_38 = total_weighte_layer_n_array[38];
          b_total_weighte_layer_n_39 = total_weighte_layer_n_array[39];
          b_total_weighte_layer_n_40 = total_weighte_layer_n_array[40];
          b_total_weighte_layer_n_41 = total_weighte_layer_n_array[41];
          b_total_weighte_layer_n_42 = total_weighte_layer_n_array[42];
          b_total_weighte_layer_n_43 = total_weighte_layer_n_array[43];
          b_total_weighte_layer_n_44 = total_weighte_layer_n_array[44];
          b_total_weighte_layer_n_45 = total_weighte_layer_n_array[45];
          b_total_weighte_layer_n_46 = total_weighte_layer_n_array[46];
          b_total_weighte_layer_n_47 = total_weighte_layer_n_array[47];
          b_total_weighte_layer_n_48 = total_weighte_layer_n_array[48];
          b_total_weighte_layer_n_49 = total_weighte_layer_n_array[49];
          b_total_weighte_layer_n_50 = total_weighte_layer_n_array[50];
          b_total_weighte_layer_n_51 = total_weighte_layer_n_array[51];
          b_total_weighte_layer_n_52 = total_weighte_layer_n_array[52];
          b_total_weighte_layer_n_53 = total_weighte_layer_n_array[53];
          b_total_weighte_layer_n_54 = total_weighte_layer_n_array[54];
          b_total_weighte_layer_n_55 = total_weighte_layer_n_array[55];

          b_total_shower_weighte_max_layer = total_weighte_ilayermax;
          b_total_shower_weighte_start_layer = total_weighte_ilayerstart;
          b_total_shower_weighte_end_layer = total_weighte_ilayerend;
          b_total_shower_weighte_start_10_layer = total_weighte_ilayerstart_10;
          b_total_shower_weighte_end_10_layer = total_weighte_ilayerend_10;
          b_total_shower_weighte_average = total_weighte_shower_averagevalue;
          b_total_shower_weighte_max = total_weighte_shower_maxvalue;

          b_total_bar_x_layer_0 = total_bar_layer_array[0][0];
          b_total_bar_x_layer_1 = total_bar_layer_array[1][0];
          b_total_bar_x_layer_2 = total_bar_layer_array[2][0];
          b_total_bar_x_layer_3 = total_bar_layer_array[3][0];
          b_total_bar_x_layer_4 = total_bar_layer_array[4][0];
          b_total_bar_x_layer_5 = total_bar_layer_array[5][0];
          b_total_bar_x_layer_6 = total_bar_layer_array[6][0];
          b_total_bar_x_layer_7 = total_bar_layer_array[7][0];
          b_total_bar_x_layer_8 = total_bar_layer_array[8][0];
          b_total_bar_x_layer_9 = total_bar_layer_array[9][0];
          b_total_bar_x_layer_10 = total_bar_layer_array[10][0];
          b_total_bar_x_layer_11 = total_bar_layer_array[11][0];
          b_total_bar_x_layer_12 = total_bar_layer_array[12][0];
          b_total_bar_x_layer_13 = total_bar_layer_array[13][0];
          b_total_bar_x_layer_14 = total_bar_layer_array[14][0];

	  b_total_bar_y_layer_0 = total_bar_layer_array[0][1];
          b_total_bar_y_layer_1 = total_bar_layer_array[1][1];
          b_total_bar_y_layer_2 = total_bar_layer_array[2][1];
          b_total_bar_y_layer_3 = total_bar_layer_array[3][1];
          b_total_bar_y_layer_4 = total_bar_layer_array[4][1];
          b_total_bar_y_layer_5 = total_bar_layer_array[5][1];
          b_total_bar_y_layer_6 = total_bar_layer_array[6][1];
          b_total_bar_y_layer_7 = total_bar_layer_array[7][1];
          b_total_bar_y_layer_8 = total_bar_layer_array[8][1];
          b_total_bar_y_layer_9 = total_bar_layer_array[9][1];
          b_total_bar_y_layer_10 = total_bar_layer_array[10][1];
          b_total_bar_y_layer_11 = total_bar_layer_array[11][1];
          b_total_bar_y_layer_12 = total_bar_layer_array[12][1];
          b_total_bar_y_layer_13 = total_bar_layer_array[13][1];
          b_total_bar_y_layer_14 = total_bar_layer_array[14][1];

          b_total_bar_r_layer_0 = total_bar_layer_array[0][2];
          b_total_bar_r_layer_1 = total_bar_layer_array[1][2];
          b_total_bar_r_layer_2 = total_bar_layer_array[2][2];
          b_total_bar_r_layer_3 = total_bar_layer_array[3][2];
          b_total_bar_r_layer_4 = total_bar_layer_array[4][2];
          b_total_bar_r_layer_5 = total_bar_layer_array[5][2];
          b_total_bar_r_layer_6 = total_bar_layer_array[6][2];
          b_total_bar_r_layer_7 = total_bar_layer_array[7][2];
          b_total_bar_r_layer_8 = total_bar_layer_array[8][2];
          b_total_bar_r_layer_9 = total_bar_layer_array[9][2];
          b_total_bar_r_layer_10 = total_bar_layer_array[10][2];
          b_total_bar_r_layer_11 = total_bar_layer_array[11][2];
          b_total_bar_r_layer_12 = total_bar_layer_array[12][2];


	  b_nhit_ratio = ratio_nhit;
	  b_sume_ratio = ratio_sume;
	  b_weighte_ratio = ratio_weighte;

	  outtree->Fill();

	  // Clear memory
	  ecal_hit_energy->clear();
	  ecal_hit_x->clear();
	  ecal_hit_y->clear();
	  ecal_hit_z->clear();
	  ecal_hit_isMasked->clear();
	  ecal_hit_slab->clear();
	  
	  hcal_hit_energy->clear();
          hcal_hit_x->clear();
          hcal_hit_y->clear();
          hcal_hit_z->clear();
          hcal_hit_isMasked->clear();
          hcal_hit_slab->clear();
	
        }
    }
    // End of writting histos	
    //Write objects
    f.WriteTObject(outtree);
        
    f.Close();
	
}
