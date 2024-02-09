#include <iostream>
#include <vector>

#include "conf_struct.h"

#define N_LAYER 15
#define N_CHIP 16
#define N_CHAN 64

void plotmap (string filename="Run_Settings/Run_Settings_90320_e-_10.0GeV.txt") {
    
  TString map_fev10="fev10_chip_channel_x_y_mapping.txt";
  TString map_cob  ="fev11_cob_rotate_chip_channel_x_y_mapping.txt";

  // Read the mapping
  map< string, vector<vector<int>> > fev10_xy = ReadMap(map_fev10);
  map< string, vector<vector<int>> > cob_xy   = ReadMap(map_cob);

  // Read the run settings
  read_configuration_file(filename,false);

  // Plot the map
  TH2F* mapxy =new TH2F("mapxy","map-xy; x; y",32,-90,90,32,-90,90);
  TH2F* mapxy_chip =new TH2F("mapxy_chip","map-xy; x; y",32,-90,90,32,-90,90);
  int fev10_layer = 4;
  for( int i_chip=0; i_chip < N_CHIP; i_chip++) {
    for( int i_chan=0; i_chan < N_CHAN; i_chan++) {

      if( detector.slab[0][fev10_layer].asu[0].skiroc[i_chip].mask[i_chan] ){
        mapxy->Fill(fev10_xy["X"].at(i_chip).at(i_chan),fev10_xy["Y"].at(i_chip).at(i_chan),i_chan);
      }
      mapxy_chip->Fill(fev10_xy["X"].at(i_chip).at(i_chan),fev10_xy["Y"].at(i_chip).at(i_chan),i_chip+1);

    }
  }

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  mapxy_chip->Draw("colz");
  mapxy->Draw("textsame");




}