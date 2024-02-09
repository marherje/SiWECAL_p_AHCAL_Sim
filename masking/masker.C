#include <iostream>
#include <vector>
#include "conf_struct.h"

#define N_LAYER 15
#define N_CHIP 16
#define N_CHAN 64

void masker (string filename="Run_Settings/Run_Settings_90320_e-_10.0GeV.txt", string in_filename="/data_ilc/flc/jimenez/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_conf6_e-_10GeV_5kevt_build.root") {

  // Read the run settings
  read_configuration_file(filename,false);

  // Apply the mask in a new branch

  TFile *oldfile = new TFile(in_filename.c_str());
  TTree *oldtree = (TTree*)oldfile->Get("ecal");

  // Activate all branches but one -> leads to a segfault!
  oldtree->SetBranchStatus("*", 0);
  auto branchnames = {"event", "spill", "cycle",
                      "bcid", "bcid_first_sca_full", "bcid_merge_end",
                      "id_run", "id_dat",
                      "nhit_slab", "nhit_chip", "nhit_chan", "nhit_len",
                      "sum_energy", "sum_energy_lg",
                      "hit_slab", "hit_chip", "hit_chan", "hit_sca",
                      "hit_adc_high", "hit_adc_low", "hit_n_scas_filled",
                      "hit_isHit", "hit_isCommissioned", "hit_energy",
                      "hit_energy_w", "hit_energy_lg",
                      "hit_x", "hit_y", "hit_z"};
  for (auto activeBranchName : branchnames) {
    oldtree->SetBranchStatus(activeBranchName, 1);
  }

  // Create a new file + a clone of old tree in new file
  size_t len = in_filename.find(".root");
  string new_filename = in_filename.substr(0, len) + "_masked.root";
  Long64_t nentries = oldtree->GetEntries();
  TFile *newfile = new TFile(new_filename.c_str(), "recreate");
  auto newtree = oldtree->CloneTree(nentries);

  // Cycle in old ttree, check and mask
  vector<int> *hit_chip = 0;
  vector<int> *hit_chan = 0;
  vector<int> *hit_slab = 0;
  vector<int> *hit_isMasked = 0;
  
  TBranch *bhit_chip = 0;
  TBranch *bhit_chan = 0;
  TBranch *bhit_slab = 0;
  TBranch *bhit_isMasked = 0;
  oldtree->SetBranchAddress("hit_chip", &hit_chip, &bhit_chip);
  oldtree->SetBranchAddress("hit_chan", &hit_chan, &bhit_chan);
  oldtree->SetBranchAddress("hit_slab", &hit_slab, &bhit_slab);

  newtree->Branch("hit_isMasked", &hit_isMasked);
  newtree->SetBranchAddress("hit_isMasked", &hit_isMasked, &bhit_isMasked);

  for (int i_event = 0; i_event < nentries; i_event++) {
    oldtree->GetEntry(i_event);
    for (int j = 0; j < hit_chip->size(); j++) {
      hit_isMasked->push_back(detector.slab[0][hit_slab->at(j)].asu[0].skiroc[hit_chip->at(j)].mask[hit_chan->at(j)]);
    }
    bhit_isMasked->Fill();
    hit_isMasked->clear();
    hit_chip->clear();
    hit_chan->clear();
    hit_slab->clear();
  }
  newfile->cd();
  newtree->Write();
  delete oldfile;
  delete newfile;

  exit(0);

}