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

#define N_ENERGIES 18

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

float moliere(vector<float> * hit_energy,
              vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z,
              vector<int> * hit_isMasked, bool masked=false, float containment = 0.95) {

    vector<float> hit_rs;
    vector<float> hit_es;
    float mol_rad = 0.;
    float sume = 0.;
    float wx = 0.; float wy = 0.; float wz = 0.;
    float r;

    for (int j = 0; j < hit_energy->size(); j++){
        if (masked && hit_isMasked->at(j) == 1) continue;
        sume += hit_energy->at(j);
        wx += hit_x->at(j) * hit_energy->at(j);
        wy += hit_y->at(j) * hit_energy->at(j);
        wz += hit_z->at(j) * hit_energy->at(j);
    }

    float bary_x = wx / sume;  
    float bary_y = wy / sume;  

    for (int j = 0; j < hit_energy->size(); j++){
        if (masked && hit_isMasked->at(j) == 1) continue;
        r = pow(pow((hit_x->at(j) - bary_x) , 2) + pow((hit_y->at(j) - bary_y), 2), 0.5);
        hit_rs.push_back(r);
    }
    for (auto k: sort_indexes(hit_rs)) hit_es.push_back(hit_energy->at(k));
    
    sort(hit_rs.begin(), hit_rs.end());

    float mol_e = 0.;
    int mol_i;
    
    for (mol_i = 0; mol_i < hit_rs.size(); mol_i++) {
        mol_e += hit_es.at(mol_i);
        if (mol_e >= containment * sume) break; 
    }
    
    return hit_rs.at(mol_i-1);      

}

void fit_res(TH1 * h, double &mu, double &sig, double &res){
    h->Fit("gaus", "q");
    mu = h->GetFunction("gaus")->GetParameter("Mean");
    sig = h->GetFunction("gaus")->GetParameter("Sigma");
    res = h->GetFunction("gaus")->GetParameter("Sigma") / h->GetFunction("gaus")->GetParameter("Mean");
    return;
}

void graph_setup_add(TGraph *g, string title, Color_t color){
    g->SetTitle(title.c_str());
    g->SetLineColor(color);
    g->SetLineWidth(3);
    return;
}

void resolutions (bool transformed = true) {
    
    // double energies[N_ENERGIES] = {1, 2, 5, 8, 10, 20, 40, 60, 80, 100, 150};
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
    for (int j = 0; j < N_ENERGIES; j++) filenames[j] = base_path + "output_LCIO2Build_TB2022-06_e-_" + to_string((int)round(energies[j])) +  "GeV.root";
    
    // for (int j = 0; j < N_ENERGIES; j++) filenames[j] = base_path + "CONF11/build/ECAL_QGSP_BERT_conf8_e-_" + to_string((int)round(energies[j])) +  "GeV_5kevt_build_masked.root";
    double test_zeros[N_ENERGIES] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // output_LCIO2Build_e_200GeV.root
    TVectorD zeros(N_ENERGIES, test_zeros);
    
    TVectorD mu_nhit(N_ENERGIES), mu_nhit_masked(N_ENERGIES), mu_sume(N_ENERGIES), mu_sume_masked(N_ENERGIES), mu_weight(N_ENERGIES), mu_weight_masked(N_ENERGIES);
    TVectorD sig_nhit(N_ENERGIES), sig_nhit_masked(N_ENERGIES), sig_sume(N_ENERGIES), sig_sume_masked(N_ENERGIES), sig_weight(N_ENERGIES), sig_weight_masked(N_ENERGIES);
    TVectorD res_nhit(N_ENERGIES), res_nhit_masked(N_ENERGIES), res_sume(N_ENERGIES), res_sume_masked(N_ENERGIES), res_weight(N_ENERGIES), res_weight_masked(N_ENERGIES);

    TVectorD mol(N_ENERGIES), mol_masked(N_ENERGIES);
    TVectorD mol_sig(N_ENERGIES), mol_sig_masked(N_ENERGIES);
    
    // Output filename
    TFile f("resolution_result.root", "recreate");

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
        
        TH1I *h_nhit = new TH1I(("NumHits_" + e_str).c_str(), ("Number of hits " + e_str).c_str(), 4000, 0, 4000);
        TH1I *h_nhit_masked = new TH1I(("NumHitsMask_" + e_str).c_str(), ("Number of hits, masked " + e_str).c_str(), 4000, 0, 4000);
        TH1F *h_sume = new TH1F(("SumEnergy_" + e_str).c_str(), ("Sum Energy " + e_str).c_str(), 500, 0, 15000);
        TH1F *h_sume_masked = new TH1F(("SumEnergyMask_" + e_str).c_str(), ("Sum Energy, masked " + e_str).c_str(), 500, 0, 15000);
        TH1F *h_weight = new TH1F(("WSumEnergy_" + e_str).c_str(), ("W Sum Energy " + e_str).c_str(), 500, 0, 70000);
        TH1F *h_weight_masked = new TH1F(("WSumEnergyMask_" + e_str).c_str(), ("W Sum Energy, masked " + e_str).c_str(), 500, 0, 70000);
        
        // Moliere histos
        TH1F *h_mol = new TH1F(("Radius95_" + e_str).c_str(), ("Radius containing 95% of energy" + e_str).c_str(), 100, 0., 100.);
        TH1F *h_mol_masked = new TH1F(("Radius95Masked_" + e_str).c_str(), ("Radius containing 95% of energy Masked " + e_str).c_str(), 100, 0., 100.);

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

	    //cout<<"debug 5"<<endl;
	    
            // Fill Moliere radii histograms
            //h_mol->Fill(moliere(hit_energy, hit_x, hit_y, hit_z, hit_isMasked));
            //h_mol_masked->Fill(moliere(hit_energy, hit_x, hit_y, hit_z, hit_isMasked, true));
            
	    //cout<<"debug 6"<<endl;
	    
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
        
	/*
        h_mol->Fit("gaus", "q");
        mol[i_energy] = h_mol->GetFunction("gaus")->GetParameter("Mean");
        mol_sig[i_energy] = h_mol->GetFunction("gaus")->GetParameter("Sigma")/nentries;
        
        h_mol_masked->Fit("gaus", "q");
        mol_masked[i_energy] = h_mol_masked->GetFunction("gaus")->GetParameter("Mean");
        mol_sig_masked[i_energy] = h_mol_masked->GetFunction("gaus")->GetParameter("Sigma")/nentries;
	*/
	
        f.WriteTObject(h_nhit);    f.WriteTObject(h_nhit_masked);
        f.WriteTObject(h_sume);    f.WriteTObject(h_sume_masked);
        f.WriteTObject(h_weight);  f.WriteTObject(h_weight_masked);
        f.WriteTObject(h_mol);     f.WriteTObject(h_mol_masked);

    }
    
    f.WriteObject(&mu_nhit, "mu_nhit");                     f.WriteObject(&sig_nhit, "sig_nhit");                   f.WriteObject(&res_nhit, "res_nhit");
    f.WriteObject(&mu_nhit_masked, "mu_nhit_masked");       f.WriteObject(&sig_nhit_masked, "sig_nhit_masked");     f.WriteObject(&res_nhit_masked, "res_nhit_masked");
    f.WriteObject(&mu_sume, "mu_sume");                     f.WriteObject(&sig_sume, "sig_sume");                   f.WriteObject(&res_sume, "res_sume");
    f.WriteObject(&mu_sume_masked, "mu_sume_masked");       f.WriteObject(&sig_sume_masked, "sig_sume_masked");     f.WriteObject(&res_sume_masked, "res_sume_masked");
    f.WriteObject(&mu_weight, "mu_weight");                 f.WriteObject(&sig_weight, "sig_weight");               f.WriteObject(&res_weight, "res_weight");
    f.WriteObject(&mu_weight_masked, "mu_weight_masked");   f.WriteObject(&sig_weight_masked, "sig_weight_masked"); f.WriteObject(&res_weight_masked, "res_weight_masked");
    f.WriteObject(&mol, "mol");   f.WriteObject(&mol_masked, "mol_masked"); f.WriteObject(&mol_sig, "mol_sig");   f.WriteObject(&mol_sig_masked, "mol_sig_masked");
    f.WriteObject(&energies, "energies");
    f.WriteObject(&energies_tr, "energies_tr");
    f.WriteObject(&W_thicknesses, "W_thicknesses");
    f.Close();

}
