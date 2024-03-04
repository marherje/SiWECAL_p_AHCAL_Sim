#define N_ENERGIES 18
#define N_ECAL_LAYERS 15
// We add the -1 layer

void addCaliceLogo(bool WIP = true){
  TImage *img = TImage::Open("style/CALICELogo_18pc.png");
  img->SetConstRatio(kTRUE);
  img->SetImageCompression(0);
  TPad *p1 = new TPad("img", "img", 0.80, 0.905, 1.0, 1.0);
  p1->Draw();
  p1->cd();
  img->Draw();

  if(WIP == true){
    TPad *p2 = new TPad("img", "img", 0.01, 0.001, 1.0, 0.2);
    p1->cd();
    p2->Draw();
    p2->cd();
    TText* t = new TText(0.15,0.3,"Work in progress");
    t->SetTextColor(kRed);
    t->SetTextSize(0.99);
    t->Draw();
  }
  
  return;
  }

void graph_setup_add(TGraph *g, string title, Color_t color){
    g->SetTitle(title.c_str());
    g->SetLineColor(color);
    g->SetLineWidth(2);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(21);
    g->SetMarkerSize(0.8);
    return;
}

TF1 *linearfit(TGraph *graph, int N, Color_t color){
  double x_axis_d[N];
  double y_axis_d[N];
  
  for(int i=0; i<N;i++){
    x_axis_d[i]=graph->GetX()[i];
    y_axis_d[i]=graph->GetY()[i];
  }
  
  TGraph *plot2fit = new TGraph(N,x_axis_d,y_axis_d);
  TF1 *f = new TF1("f", "[1]*x +[0]",0.05,0.75);
  f->SetLineColor(color);
  plot2fit->Fit(f);
  return f;
}

TF1 *linearfitErrors(TGraphErrors *grapherrors, int N, Color_t color){
  double x_axis_d[N];
  double y_axis_d[N];

  for(int i=0; i<N;i++){
    x_axis_d[i]=grapherrors->GetX()[i];
    y_axis_d[i]=grapherrors->GetY()[i];
    cout<<"(x,y): ("<<x_axis_d[i]<<","<<y_axis_d[i]<<")"<<endl;
  }

  TGraph *plot2fit = new TGraph(N,x_axis_d,y_axis_d);
  TF1 *f = new TF1("f", "[1]*x +[0]",1,210);
  f->SetLineColor(color);
  plot2fit->Fit(f);
  return f;
}

void res_mol_plots(string particle, bool transformed = true, bool save = false){
    
    TString filename = "results_folder/resolution_"+particle+"_result.root";
    TFile *file = new TFile(filename, "read");
    // TVectorD * energies = file->GetObject("energies");
    TString savetrans;
    if(transformed == false) savetrans="Linear_"+particle+"_";
    else savetrans="InvSq_"+particle+"_";
    
    TVectorD *energies;
    TVectorD *energies_tr;
    
    TVectorD *res_nhit, *mu_nhit, *sig_nhit;
    TVectorD *res_nhit_masked, *mu_nhit_masked, *sig_nhit_masked;
    TVectorD *res_sume, *mu_sume, *sig_sume;
    TVectorD *res_sume_masked, *mu_sume_masked, *sig_sume_masked;
    TVectorD *res_weight, *mu_weight, *sig_weight;
    TVectorD *res_weight_masked, *mu_weight_masked, *sig_weight_masked;
    
    TVectorD *mu_mol, *mu_mol_masked, *sig_mol, *sig_mol_masked;

    TVectorD *mu_nhit_layer_0, *sig_nhit_layer_0;
    TVectorD *mu_nhit_layer_1, *sig_nhit_layer_1;
    TVectorD *mu_nhit_layer_2, *sig_nhit_layer_2;
    TVectorD *mu_nhit_layer_3, *sig_nhit_layer_3;
    TVectorD *mu_nhit_layer_4, *sig_nhit_layer_4;
    TVectorD *mu_nhit_layer_5, *sig_nhit_layer_5;
    TVectorD *mu_nhit_layer_6, *sig_nhit_layer_6;
    TVectorD *mu_nhit_layer_7, *sig_nhit_layer_7;
    TVectorD *mu_nhit_layer_8, *sig_nhit_layer_8;
    TVectorD *mu_nhit_layer_9, *sig_nhit_layer_9;
    TVectorD *mu_nhit_layer_10, *sig_nhit_layer_10;
    TVectorD *mu_nhit_layer_11, *sig_nhit_layer_11;
    TVectorD *mu_nhit_layer_12, *sig_nhit_layer_12;
    TVectorD *mu_nhit_layer_13, *sig_nhit_layer_13;
    TVectorD *mu_nhit_layer_14, *sig_nhit_layer_14;

    TVectorD *shower_nhit_max_layer, *shower_nhit_start_layer, *shower_nhit_end_layer, *shower_nhit_max, *shower_nhit_average;

    TVectorD *mu_nhit_layer_w_0, *sig_nhit_layer_w_0;
    TVectorD *mu_nhit_layer_w_1, *sig_nhit_layer_w_1;
    TVectorD *mu_nhit_layer_w_2, *sig_nhit_layer_w_2;
    TVectorD *mu_nhit_layer_w_3, *sig_nhit_layer_w_3;
    TVectorD *mu_nhit_layer_w_4, *sig_nhit_layer_w_4;
    TVectorD *mu_nhit_layer_w_5, *sig_nhit_layer_w_5;
    TVectorD *mu_nhit_layer_w_6, *sig_nhit_layer_w_6;
    TVectorD *mu_nhit_layer_w_7, *sig_nhit_layer_w_7;
    TVectorD *mu_nhit_layer_w_8, *sig_nhit_layer_w_8;
    TVectorD *mu_nhit_layer_w_9, *sig_nhit_layer_w_9;
    TVectorD *mu_nhit_layer_w_10, *sig_nhit_layer_w_10;
    TVectorD *mu_nhit_layer_w_11, *sig_nhit_layer_w_11;
    TVectorD *mu_nhit_layer_w_12, *sig_nhit_layer_w_12;
    TVectorD *mu_nhit_layer_w_13, *sig_nhit_layer_w_13;
    TVectorD *mu_nhit_layer_w_14, *sig_nhit_layer_w_14;


    // New variables
    // Index    0 1 2 3  4  5  6  7  8  9 10 11 12  13  14  15  16  17
    // Energies 2 4 6 8 10 20 30 40 50 60 70 80 90 100 125 150 175 200
    double showervalues[N_ENERGIES][N_ECAL_LAYERS];
    double showervalues_w[N_ENERGIES][N_ECAL_LAYERS];

    TVectorD *layeraxis = new TVectorD(8);
    double layeraxisvalues[N_ECAL_LAYERS] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.};
    TVectorD layeraxisref(N_ECAL_LAYERS, layeraxisvalues);
    layeraxis = &layeraxisref;
    
    TVectorD *shower_profile_2 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_2[N_ECAL_LAYERS];
    TVectorD *shower_profile_4 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_4[N_ECAL_LAYERS];
    TVectorD *shower_profile_6 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_6[N_ECAL_LAYERS];
    TVectorD *shower_profile_10 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_10[N_ECAL_LAYERS];
    TVectorD *shower_profile_20 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_20[N_ECAL_LAYERS];
    TVectorD *shower_profile_40 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_40[N_ECAL_LAYERS];
    TVectorD *shower_profile_80 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_80[N_ECAL_LAYERS];
    TVectorD *shower_profile_150 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_150[N_ECAL_LAYERS];

    TVectorD *shower_profile_w_2 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_w_2[N_ECAL_LAYERS];
    TVectorD *shower_profile_w_4 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_w_4[N_ECAL_LAYERS];
    TVectorD *shower_profile_w_6 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_w_6[N_ECAL_LAYERS];
    TVectorD *shower_profile_w_10 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_w_10[N_ECAL_LAYERS];
    TVectorD *shower_profile_w_20 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_w_20[N_ECAL_LAYERS];
    TVectorD *shower_profile_w_40 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_w_40[N_ECAL_LAYERS];
    TVectorD *shower_profile_w_80 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_w_80[N_ECAL_LAYERS];
    TVectorD *shower_profile_w_150 = new TVectorD(N_ECAL_LAYERS);
    double showervalues_w_150[N_ECAL_LAYERS];

    TVectorD *zeros = new TVectorD(N_ENERGIES);
    double zerovalues[N_ENERGIES] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}; 
    TVectorD zerovectorref(N_ENERGIES, zerovalues);
    zeros = &zerovectorref;

    // +1 fr the -1 layer
    TVectorD *zeroslayer = new TVectorD(N_ECAL_LAYERS);
    double zerovalueslayer[N_ECAL_LAYERS] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    TVectorD zerovectorreflayer(N_ECAL_LAYERS, zerovalueslayer);
    zeroslayer = &zerovectorreflayer;
    
    // Load objects
    file->GetObject("energies", energies);
    file->GetObject("energies_tr", energies_tr);
    
    file->GetObject("res_nhit", res_nhit);
    file->GetObject("res_nhit_masked", res_nhit_masked);
    file->GetObject("res_sume", res_sume);
    file->GetObject("res_sume_masked", res_sume_masked);
    file->GetObject("res_weight", res_weight);
    file->GetObject("res_weight_masked", res_weight_masked);
    
    file->GetObject("mu_nhit", mu_nhit);
    file->GetObject("mu_nhit_masked", mu_nhit_masked);
    file->GetObject("mu_sume", mu_sume);
    file->GetObject("mu_sume_masked", mu_sume_masked);
    file->GetObject("mu_weight", mu_weight);
    file->GetObject("mu_weight_masked", mu_weight_masked);
    
    file->GetObject("sig_nhit", sig_nhit);
    file->GetObject("sig_nhit_masked", sig_nhit_masked);
    file->GetObject("sig_sume", sig_sume);
    file->GetObject("sig_sume_masked", sig_sume_masked);
    file->GetObject("sig_weight", sig_weight);
    file->GetObject("sig_weight_masked", sig_weight_masked);
    
    file->GetObject("mu_mol", mu_mol);
    file->GetObject("mu_mol_masked", mu_mol_masked);
    file->GetObject("sig_mol", sig_mol);
    file->GetObject("sig_mol_masked", sig_mol_masked);

    file->GetObject("mu_nhit_layer_0", mu_nhit_layer_0);
    file->GetObject("mu_nhit_layer_1", mu_nhit_layer_1);
    file->GetObject("mu_nhit_layer_2", mu_nhit_layer_2);
    file->GetObject("mu_nhit_layer_3", mu_nhit_layer_3);
    file->GetObject("mu_nhit_layer_4", mu_nhit_layer_4);
    file->GetObject("mu_nhit_layer_5", mu_nhit_layer_5);
    file->GetObject("mu_nhit_layer_6", mu_nhit_layer_6);
    file->GetObject("mu_nhit_layer_7", mu_nhit_layer_7);
    file->GetObject("mu_nhit_layer_8", mu_nhit_layer_8);
    file->GetObject("mu_nhit_layer_9", mu_nhit_layer_9);
    file->GetObject("mu_nhit_layer_10", mu_nhit_layer_10);
    file->GetObject("mu_nhit_layer_11", mu_nhit_layer_11);
    file->GetObject("mu_nhit_layer_12", mu_nhit_layer_12);
    file->GetObject("mu_nhit_layer_13", mu_nhit_layer_13);
    file->GetObject("mu_nhit_layer_14", mu_nhit_layer_14);
    file->GetObject("sig_nhit_layer_0", sig_nhit_layer_0);
    file->GetObject("sig_nhit_layer_1", sig_nhit_layer_1);
    file->GetObject("sig_nhit_layer_2", sig_nhit_layer_2);
    file->GetObject("sig_nhit_layer_3", sig_nhit_layer_3);
    file->GetObject("sig_nhit_layer_4", sig_nhit_layer_4);
    file->GetObject("sig_nhit_layer_5", sig_nhit_layer_5);
    file->GetObject("sig_nhit_layer_6", sig_nhit_layer_6);
    file->GetObject("sig_nhit_layer_7", sig_nhit_layer_7);
    file->GetObject("sig_nhit_layer_8", sig_nhit_layer_8);
    file->GetObject("sig_nhit_layer_9", sig_nhit_layer_9);
    file->GetObject("sig_nhit_layer_10", sig_nhit_layer_10);
    file->GetObject("sig_nhit_layer_11", sig_nhit_layer_11);
    file->GetObject("sig_nhit_layer_12", sig_nhit_layer_12);
    file->GetObject("sig_nhit_layer_13", sig_nhit_layer_13);
    file->GetObject("sig_nhit_layer_14", sig_nhit_layer_14);

    file->GetObject("mu_nhit_layer_w_0", mu_nhit_layer_w_0);
    file->GetObject("mu_nhit_layer_w_1", mu_nhit_layer_w_1);
    file->GetObject("mu_nhit_layer_w_2", mu_nhit_layer_w_2);
    file->GetObject("mu_nhit_layer_w_3", mu_nhit_layer_w_3);
    file->GetObject("mu_nhit_layer_w_4", mu_nhit_layer_w_4);
    file->GetObject("mu_nhit_layer_w_5", mu_nhit_layer_w_5);
    file->GetObject("mu_nhit_layer_w_6", mu_nhit_layer_w_6);
    file->GetObject("mu_nhit_layer_w_7", mu_nhit_layer_w_7);
    file->GetObject("mu_nhit_layer_w_8", mu_nhit_layer_w_8);
    file->GetObject("mu_nhit_layer_w_9", mu_nhit_layer_w_9);
    file->GetObject("mu_nhit_layer_w_10", mu_nhit_layer_w_10);
    file->GetObject("mu_nhit_layer_w_11", mu_nhit_layer_w_11);
    file->GetObject("mu_nhit_layer_w_12", mu_nhit_layer_w_12);
    file->GetObject("mu_nhit_layer_w_13", mu_nhit_layer_w_13);
    file->GetObject("mu_nhit_layer_w_14", mu_nhit_layer_w_14);
    file->GetObject("sig_nhit_layer_w_0", sig_nhit_layer_w_0);
    file->GetObject("sig_nhit_layer_w_1", sig_nhit_layer_w_1);
    file->GetObject("sig_nhit_layer_w_2", sig_nhit_layer_w_2);
    file->GetObject("sig_nhit_layer_w_3", sig_nhit_layer_w_3);
    file->GetObject("sig_nhit_layer_w_4", sig_nhit_layer_w_4);
    file->GetObject("sig_nhit_layer_w_5", sig_nhit_layer_w_5);
    file->GetObject("sig_nhit_layer_w_6", sig_nhit_layer_w_6);
    file->GetObject("sig_nhit_layer_w_7", sig_nhit_layer_w_7);
    file->GetObject("sig_nhit_layer_w_8", sig_nhit_layer_w_8);
    file->GetObject("sig_nhit_layer_w_9", sig_nhit_layer_w_9);
    file->GetObject("sig_nhit_layer_w_10", sig_nhit_layer_w_10);
    file->GetObject("sig_nhit_layer_w_11", sig_nhit_layer_w_11);
    file->GetObject("sig_nhit_layer_w_12", sig_nhit_layer_w_12);
    file->GetObject("sig_nhit_layer_w_13", sig_nhit_layer_w_13);
    file->GetObject("sig_nhit_layer_w_14", sig_nhit_layer_w_14);

    file->GetObject("shower_nhit_max_layer", shower_nhit_max_layer);
    file->GetObject("shower_nhit_start_layer", shower_nhit_start_layer);
    file->GetObject("shower_nhit_end_layer", shower_nhit_end_layer);
    file->GetObject("shower_nhit_max", shower_nhit_max);
    file->GetObject("shower_nhit_average", shower_nhit_average);

    // Load all the showers
    for(int ienergy=0; ienergy < N_ENERGIES; ienergy++){
      showervalues[ienergy][0]=((*mu_nhit_layer_0))[ienergy];
      showervalues[ienergy][1]=((*mu_nhit_layer_1))[ienergy];
      showervalues[ienergy][2]=((*mu_nhit_layer_2))[ienergy];
      showervalues[ienergy][3]=((*mu_nhit_layer_3))[ienergy];
      showervalues[ienergy][4]=((*mu_nhit_layer_4))[ienergy];
      showervalues[ienergy][5]=((*mu_nhit_layer_5))[ienergy];
      showervalues[ienergy][6]=((*mu_nhit_layer_6))[ienergy];
      showervalues[ienergy][7]=((*mu_nhit_layer_7))[ienergy];
      showervalues[ienergy][8]=((*mu_nhit_layer_8))[ienergy];
      showervalues[ienergy][9]=((*mu_nhit_layer_9))[ienergy];
      showervalues[ienergy][10]=((*mu_nhit_layer_10))[ienergy];
      showervalues[ienergy][11]=((*mu_nhit_layer_11))[ienergy];
      showervalues[ienergy][12]=((*mu_nhit_layer_12))[ienergy];
      showervalues[ienergy][13]=((*mu_nhit_layer_13))[ienergy];
      showervalues[ienergy][14]=((*mu_nhit_layer_14))[ienergy];
      
      showervalues_w[ienergy][0]=((*mu_nhit_layer_w_0))[ienergy];
      showervalues_w[ienergy][1]=((*mu_nhit_layer_w_1))[ienergy];
      showervalues_w[ienergy][2]=((*mu_nhit_layer_w_2))[ienergy];
      showervalues_w[ienergy][3]=((*mu_nhit_layer_w_3))[ienergy];
      showervalues_w[ienergy][4]=((*mu_nhit_layer_w_4))[ienergy];
      showervalues_w[ienergy][5]=((*mu_nhit_layer_w_5))[ienergy];
      showervalues_w[ienergy][6]=((*mu_nhit_layer_w_6))[ienergy];
      showervalues_w[ienergy][7]=((*mu_nhit_layer_w_7))[ienergy];
      showervalues_w[ienergy][8]=((*mu_nhit_layer_w_8))[ienergy];
      showervalues_w[ienergy][9]=((*mu_nhit_layer_w_9))[ienergy];
      showervalues_w[ienergy][10]=((*mu_nhit_layer_w_10))[ienergy];
      showervalues_w[ienergy][11]=((*mu_nhit_layer_w_11))[ienergy];
      showervalues_w[ienergy][12]=((*mu_nhit_layer_w_12))[ienergy];
      showervalues_w[ienergy][13]=((*mu_nhit_layer_w_13))[ienergy];
      showervalues_w[ienergy][14]=((*mu_nhit_layer_w_14))[ienergy];
    }

    for(int ilayer=0; ilayer < N_ECAL_LAYERS; ilayer++){
      // Index    0 1 2 3  4  5  6  7  8  9 10 11 12  13  14  15  16  17                                                         
      // Energies 2 4 6 8 10 20 30 40 50 60 70 80 90 100 125 150 175 200
      showervalues_2[ilayer] = showervalues[0][ilayer];
      showervalues_4[ilayer] = showervalues[1][ilayer];
      showervalues_6[ilayer] = showervalues[2][ilayer];
      showervalues_10[ilayer] = showervalues[4][ilayer];
      showervalues_20[ilayer] = showervalues[5][ilayer];
      showervalues_40[ilayer] = showervalues[7][ilayer];
      showervalues_80[ilayer] = showervalues[11][ilayer];
      showervalues_150[ilayer] = showervalues[15][ilayer];
    
      showervalues_w_2[ilayer] = showervalues_w[0][ilayer];
      showervalues_w_4[ilayer] = showervalues_w[1][ilayer];
      showervalues_w_6[ilayer] = showervalues_w[2][ilayer];
      showervalues_w_10[ilayer] = showervalues_w[4][ilayer];
      showervalues_w_20[ilayer] = showervalues_w[5][ilayer];
      showervalues_w_40[ilayer] = showervalues_w[7][ilayer];
      showervalues_w_80[ilayer] = showervalues_w[11][ilayer];
      showervalues_w_150[ilayer] = showervalues_w[15][ilayer];
    }
    
    TVectorD showerprofileref_2(N_ECAL_LAYERS,showervalues_2);
    shower_profile_2 = &showerprofileref_2;
    TVectorD showerprofileref_4(N_ECAL_LAYERS,showervalues_4);
    shower_profile_4 = &showerprofileref_4;
    TVectorD showerprofileref_6(N_ECAL_LAYERS,showervalues_6);
    shower_profile_6 = &showerprofileref_6;
    TVectorD showerprofileref_10(N_ECAL_LAYERS,showervalues_10);
    shower_profile_10 = &showerprofileref_10;
    TVectorD showerprofileref_20(N_ECAL_LAYERS,showervalues_20);
    shower_profile_20 = &showerprofileref_20;
    TVectorD showerprofileref_40(N_ECAL_LAYERS,showervalues_40);
    shower_profile_40 = &showerprofileref_40;
    TVectorD showerprofileref_80(N_ECAL_LAYERS,showervalues_80);
    shower_profile_80 = &showerprofileref_80;
    TVectorD showerprofileref_150(N_ECAL_LAYERS,showervalues_150);
    shower_profile_150 = &showerprofileref_150;

    TVectorD showerprofileref_w_2(N_ECAL_LAYERS,showervalues_w_2);
    shower_profile_w_2 = &showerprofileref_w_2;
    TVectorD showerprofileref_w_4(N_ECAL_LAYERS,showervalues_w_4);
    shower_profile_w_4 = &showerprofileref_w_4;
    TVectorD showerprofileref_w_6(N_ECAL_LAYERS,showervalues_w_6);
    shower_profile_w_6 = &showerprofileref_w_6;
    TVectorD showerprofileref_w_10(N_ECAL_LAYERS,showervalues_w_10);
    shower_profile_w_10 = &showerprofileref_w_10;
    TVectorD showerprofileref_w_20(N_ECAL_LAYERS,showervalues_w_20);
    shower_profile_w_20 = &showerprofileref_w_20;
    TVectorD showerprofileref_w_40(N_ECAL_LAYERS,showervalues_w_40);
    shower_profile_w_40 = &showerprofileref_w_40;
    TVectorD showerprofileref_w_80(N_ECAL_LAYERS,showervalues_w_80);
    shower_profile_w_80 = &showerprofileref_w_80;
    TVectorD showerprofileref_w_150(N_ECAL_LAYERS,showervalues_w_150);
    shower_profile_w_150 = &showerprofileref_w_150;



    // Prepare energy or 1/sqrt(E)
    TVectorD *energyaxis = new TVectorD(N_ENERGIES);
    double energyvalues[N_ENERGIES];
    for(int i=0; i < N_ENERGIES; i++){
      if (transformed == true) energyvalues[i]=((*energies_tr))[i];
      else energyvalues[i]=((*energies))[i];
    }
    TVectorD energyaxisref(N_ENERGIES, energyvalues);
    energyaxis = &energyaxisref;

    
    // // string e_str = to_string((int)round(energies[i_energy])) + "GeV";
    // // Resolution
    
    TGraph *g_res_nhit = new TGraph(*energyaxis, *res_nhit);                 
    graph_setup_add(g_res_nhit, "N hits", kBlack);
    for (int i=0;i<g_res_nhit->GetN();i++) g_res_nhit->GetY()[i] *= 100;
    TF1 *fit_res_nhit = linearfit(g_res_nhit,18,kBlack);

    TGraph *g_res_nhit_masked = new TGraph(*energyaxis, *res_nhit_masked);   
    graph_setup_add(g_res_nhit_masked, "N hits (masked)", kBlue);
    for (int i=0;i<g_res_nhit_masked->GetN();i++) g_res_nhit_masked->GetY()[i] *= 100;
    
    
    TGraph *g_res_sume = new TGraph(*energyaxis, *res_sume);                 
    graph_setup_add(g_res_sume, "Sum energy", kRed);
    for (int i=0;i<g_res_sume->GetN();i++) g_res_sume->GetY()[i] *= 100;
    TF1 *fit_res_sume = linearfit(g_res_sume,18,kRed);
    
    TGraph *g_res_sume_masked = new TGraph(*energyaxis, *res_sume_masked);   
    graph_setup_add(g_res_sume_masked, "Sum energy (masked)", kGreen);
    for (int i=0;i<g_res_sume_masked->GetN();i++) g_res_sume_masked->GetY()[i] *= 100;
    
    TGraph *g_res_weight = new TGraph(*energyaxis, *res_weight);                 
    graph_setup_add(g_res_weight, "Weighted sum energy", kViolet);
    for (int i=0;i<g_res_weight->GetN();i++) g_res_weight->GetY()[i] *= 100;
    TF1 *fit_res_weight = linearfit(g_res_weight,18,kViolet);
    
    TGraph *g_res_weight_masked = new TGraph(*energyaxis, *res_weight_masked);   
    graph_setup_add(g_res_weight_masked, "Weighted sum energy (masked)", kCyan);
    for (int i=0;i<g_res_weight_masked->GetN();i++) g_res_weight_masked->GetY()[i] *= 100;
    
    // Debugging tvectord access 
    for (int j = 0; j < N_ENERGIES; j++) {cout << "This energy: " << ((*energies))[j] << " GeV, 1/sqrt(E):" << ((*energies_tr))[j] << endl;}


    //Resolution plot:
    TMultiGraph *mg = new TMultiGraph();
    
    mg->Add(g_res_nhit);
    // mg->Add(g_res_nhit_masked);
    mg->Add(g_res_sume);
    // mg->Add(g_res_sume_masked);
    mg->Add(g_res_weight);
    // mg->Add(g_res_weight_masked);
    
    auto c = new TCanvas("c", "c", 800, 800);
    mg->Draw("AP");
    mg->SetTitle("Resolution vs energy");
    if (transformed) mg->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg->GetXaxis()->SetTitle("E[GeV]");
    mg->GetYaxis()->SetTitle("#sigma(E_{meas.})/#mu(E_{meas.}) (\%)");
    if( transformed == true){
      fit_res_nhit->Draw("SAME");
      fit_res_sume->Draw("SAME");
      fit_res_weight->Draw("SAME");
    }
    TLegend *leg_r;
    if(transformed == true) leg_r= new TLegend(0.45,0.15,0.65,0.35);
    else leg_r= new TLegend(0.45,0.65,0.65,0.85);
    leg_r->SetTextSize(0.035);
    leg_r->SetTextFont(42);
    leg_r->AddEntry(g_res_nhit,"N Hits","lp");
    leg_r->AddEntry(g_res_sume,"Summed energy","lp");
    leg_r->AddEntry(g_res_weight,"Summed energy weighted","lp");
    leg_r->SetFillColor(0);
    leg_r->SetLineColor(0);
    leg_r->SetShadowColor(0);
    leg_r->Draw();
    addCaliceLogo();
    if(save == true){
    c->SaveAs("resolution_ecal_sim_"+savetrans+".eps");
    //c->SaveAs("resolution_ecal_sim_"+savetrans+".root");
    }
    
    // // Linearity plots
    // We could add sim, sim_masked and data alltogether
        
    TMultiGraph *mg_lin_nhit = new TMultiGraph();
    TMultiGraph *mg_lin_sume = new TMultiGraph();
    TMultiGraph *mg_lin_weight = new TMultiGraph();
    
    TGraphErrors *g_lin_nhit = new TGraphErrors(*energyaxis, *mu_nhit, *zeros, *sig_nhit);                 
    graph_setup_add(g_lin_nhit, "n hits", kBlack);
    TF1 *fit_lin_nhit = linearfitErrors(g_lin_nhit,18,kBlack);
    
    TGraphErrors *g_lin_nhit_masked = new TGraphErrors(*energyaxis, *mu_nhit_masked, *zeros, *sig_nhit_masked);   
    graph_setup_add(g_lin_nhit_masked, "n hits (masked)", kBlue);
    
    TGraphErrors *g_lin_sume = new TGraphErrors(*energyaxis, *mu_sume, *zeros, *sig_sume);                 
    graph_setup_add(g_lin_sume, "sum energy", kRed);
    TF1 *fit_lin_sume = linearfitErrors(g_lin_sume,18,kRed);
    
    TGraphErrors *g_lin_sume_masked = new TGraphErrors(*energyaxis, *mu_sume_masked, *zeros, *sig_sume_masked);   
    graph_setup_add(g_lin_sume_masked, "sum energy (masked)", kGreen);
    
    TGraphErrors *g_lin_weight = new TGraphErrors(*energyaxis, *mu_weight, *zeros, *sig_weight);                 
    graph_setup_add(g_lin_weight, "weighted sum energy", kViolet);
    TF1 *fit_lin_weight = linearfitErrors(g_lin_weight,18,kViolet);

    TGraphErrors *g_lin_weight_masked = new TGraphErrors(*energyaxis, *mu_weight_masked, *zeros, *sig_weight_masked);   
    graph_setup_add(g_lin_weight_masked, "weighted sum energy (weighted)", kCyan);
    

    mg_lin_nhit->Add(g_lin_nhit);
    // mg_lin_nhit->Add(g_lin_nhit_masked);
    mg_lin_sume->Add(g_lin_sume);
    // mg_lin_sume->Add(g_lin_sume_masked);
    mg_lin_weight->Add(g_lin_weight);
    // mg_lin_weight->Add(g_lin_weight_masked);
    
    auto c_lin_nhit = new TCanvas("c_lin_nhit", "c_lin_nhit", 800, 800);
    mg_lin_nhit->Draw("AP");
    mg_lin_nhit->SetTitle("Linearity (N hits)");
    if (transformed) mg_lin_nhit->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_lin_nhit->GetXaxis()->SetTitle("E[GeV]");
    mg_lin_nhit->GetYaxis()->SetTitle("Number of hits");
    mg_lin_nhit->GetYaxis()->SetMaxDigits(3);
    if( transformed == false) fit_lin_nhit->Draw("SAME");
    TLegend *leg_ln;
    if(transformed == true) leg_ln= new TLegend(0.4,0.6,0.6,0.6);
    else leg_ln= new TLegend(0.45,0.15,0.65,0.35);
    leg_ln->SetTextSize(0.035);
    leg_ln->SetTextFont(42);
    leg_ln->AddEntry(g_lin_nhit,"N Hits","lp");
    leg_ln->SetFillColor(0);
    leg_ln->SetLineColor(0);
    leg_ln->SetShadowColor(0);
    leg_ln->Draw();
    addCaliceLogo();
    if(save == true){
      c_lin_nhit->SaveAs("lin_nhits_ecal_sim_"+savetrans+".eps");
      //c_lin_nhit->SaveAs("lin_nhits_ecal_sim_"+savetrans+".root");
    }

    auto c_lin_sume = new TCanvas("c_lin_sume", "c_lin_sume", 800, 800);
    mg_lin_sume->Draw("AP");
    mg_lin_sume->SetTitle("Linearity (energy)");
    if (transformed) mg_lin_sume->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_lin_sume->GetXaxis()->SetTitle("E[GeV]");
    mg_lin_sume->GetYaxis()->SetTitle("Summed energy (MIPs)");
    mg_lin_sume->GetYaxis()->SetMaxDigits(3);
    if( transformed == false) fit_lin_sume->Draw("SAME");
    TLegend *leg_s;
    if(transformed == true) leg_s= new TLegend(0.4,0.6,0.6,0.8);
    else leg_s= new TLegend(0.45,0.15,0.65,0.35);
    leg_s->SetTextSize(0.035);
    leg_s->SetTextFont(42);
    leg_s->AddEntry(g_lin_sume,"Summed energy","lp");
    leg_s->SetFillColor(0);
    leg_s->SetLineColor(0);
    leg_s->SetShadowColor(0);
    leg_s->Draw();
    addCaliceLogo();
    if(save == true){
      c_lin_sume->SaveAs("lin_sume_ecal_sim_"+savetrans+".eps");
      //c_lin_sume->SaveAs("lin_sume_ecal_sim_"+savetrans+".root");
    }

    auto c_lin_weight = new TCanvas("c_lin_weight", "c_lin_weight", 800, 800);
    mg_lin_weight->Draw("AP");
    mg_lin_weight->SetTitle("Linearity (weighted energy)");
    if (transformed) mg_lin_weight->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_lin_weight->GetXaxis()->SetTitle("E[GeV]");
    mg_lin_weight->GetYaxis()->SetTitle("Weighted summed energy (MIPs)");
    mg_lin_weight->GetYaxis()->SetMaxDigits(3);
    if( transformed == false) fit_lin_weight->Draw("SAME");
    TLegend *leg_w;
    if(transformed == true) leg_w= new TLegend(0.4,0.6,0.6,0.8);
    else leg_w= new TLegend(0.45,0.15,0.65,0.35);
    leg_w->SetTextSize(0.035);
    leg_w->SetTextFont(42);
    leg_w->AddEntry(g_lin_weight,"Weighted summed energy","lp");
    leg_w->SetFillColor(0);
    leg_w->SetLineColor(0);
    leg_w->SetShadowColor(0);
    leg_w->Draw();
    addCaliceLogo();
    if(save == true){
      c_lin_weight->SaveAs("lin_weight_ecal_sim_"+savetrans+".eps");
      //c_lin_weight->SaveAs("lin_weight_ecal_sim_"+savetrans+".root");
    }

    
    
    // Moliere
    // TGraphErrors *g_mol = new TGraphErrors(N_ENERGIES, energies_tr, mol, zeros, mol_sig);   
    TGraphErrors *g_mol = new TGraphErrors(*energyaxis, *mu_mol, *zeros, *sig_mol);   
    graph_setup_add(g_mol, "Moliere Radius", kBlue);

    // TGraphErrors *g_mol_masked = new TGraphErrors(N_ENERGIES, energies_tr, mol_masked, zeros, mol_sig_masked);   
    TGraphErrors *g_mol_masked = new TGraphErrors(*energyaxis, *mu_mol_masked, *zeros, *sig_mol_masked);   
    graph_setup_add(g_mol_masked, "Moliere Radius (masked)", kRed);

    TMultiGraph *mg_mol = new TMultiGraph();
    mg_mol->Add(g_mol);
    // mg_mol->Add(g_mol_masked);

    auto c_mol = new TCanvas("c_mol", "c_mol", 800, 800);
    mg_mol->Draw("AP");
    mg_mol->SetTitle("Moliere radius vs energy");
    // mg_mol->GetXaxis()->SetTitle("Energy [GeV]");
    if (transformed) mg_mol->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_mol->GetXaxis()->SetTitle("E[GeV]");
    mg_mol->GetYaxis()->SetTitle("Radius (90%) [mm]");
    mg_mol->GetYaxis()->SetRangeUser(0,45);
				 
    TLegend *leg_m;
    if(transformed == true) leg_m= new TLegend(0.45,0.65,0.65,0.85);
    else leg_m= new TLegend(0.45,0.65,0.65,0.85);
    leg_m->SetTextSize(0.035);
    leg_m->SetTextFont(42);
    leg_m->AddEntry(g_mol,"Moliere radius","ap");
    leg_m->SetFillColor(0);
    leg_m->SetLineColor(0);
    leg_m->SetShadowColor(0);
    leg_m->Draw();
    addCaliceLogo();
    if(save == true){
      c_mol->SaveAs("mol_"+savetrans+".eps");
      //c_lin_weight->SaveAs("lin_weight_ecal_sim_"+savetrans+".root");                                                                                       
    }

    // Shower profile
    // Example layer 3
    TGraphErrors *g_nhits_layer_3 = new TGraphErrors(*energyaxis, *mu_nhit_layer_3, *zeros, *sig_nhit_layer_3);
    graph_setup_add(g_nhits_layer_3, "Hits in layer 3", kBlue);

    TMultiGraph *mg_nhits_layer_3 = new TMultiGraph();
    mg_nhits_layer_3->Add(g_nhits_layer_3);

    auto c_layer_3 = new TCanvas("c_layer_3", "c_layer_3", 800, 800);
    mg_nhits_layer_3->Draw("AP");
    mg_nhits_layer_3->SetTitle("Hits in layer 3");

    if (transformed) mg_nhits_layer_3->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_nhits_layer_3->GetXaxis()->SetTitle("E[GeV]");
    mg_nhits_layer_3->GetYaxis()->SetTitle("Hits");
    mg_nhits_layer_3->GetYaxis()->SetRangeUser(0,200);

    TLegend *leg_layer_3;
    if(transformed == true) leg_layer_3= new TLegend(0.45,0.65,0.65,0.85);
    else leg_layer_3= new TLegend(0.45,0.65,0.65,0.85);
    leg_layer_3->SetTextSize(0.035);
    leg_layer_3->SetTextFont(42);
    leg_layer_3->AddEntry(g_mol,"Hits in layer 3","ap");
    leg_layer_3->SetFillColor(0);
    leg_layer_3->SetLineColor(0);
    leg_layer_3->SetShadowColor(0);
    leg_layer_3->Draw();
    addCaliceLogo();
    if(save == true){
      c_layer_3->SaveAs("hits_layer_3_"+savetrans+".eps");
    }

    // Shower different energies
    TGraphErrors *g_shower_2 = new TGraphErrors(*layeraxis, *shower_profile_2 , *zeros, *zeros);
    graph_setup_add(g_shower_2, "Shower profile (2 GeV)", kBlack);
    TGraphErrors *g_shower_4 = new TGraphErrors(*layeraxis, *shower_profile_4 , *zeros, *zeros);
    graph_setup_add(g_shower_4, "Shower profile (4 GeV)", kBlue+1);
    TGraphErrors *g_shower_6 = new TGraphErrors(*layeraxis, *shower_profile_6 , *zeros, *zeros);
    graph_setup_add(g_shower_6, "Shower profile (6 GeV)", kBlue);
    TGraphErrors *g_shower_10 = new TGraphErrors(*layeraxis, *shower_profile_10 , *zeros, *zeros);
    graph_setup_add(g_shower_10, "Shower profile (10 GeV)", kAzure-4);
    TGraphErrors *g_shower_20 = new TGraphErrors(*layeraxis, *shower_profile_20 , *zeros, *zeros);
    graph_setup_add(g_shower_20, "Shower profile (20 GeV)", kAzure+8);
    TGraphErrors *g_shower_40 = new TGraphErrors(*layeraxis, *shower_profile_40 , *zeros, *zeros);
    graph_setup_add(g_shower_40, "Shower profile (40 GeV)", kRed-7);
    TGraphErrors *g_shower_80 = new TGraphErrors(*layeraxis, *shower_profile_80 , *zeros, *zeros);
    graph_setup_add(g_shower_80, "Shower profile (80 GeV)", kRed);
    TGraphErrors *g_shower_150 = new TGraphErrors(*layeraxis, *shower_profile_150 , *zeros, *zeros);
    graph_setup_add(g_shower_150, "Shower profile (150 GeV)", kRed-3);

    TMultiGraph *mg_shower_energies = new TMultiGraph();
    mg_shower_energies->Add(g_shower_2);
    mg_shower_energies->Add(g_shower_4);
    mg_shower_energies->Add(g_shower_6);
    mg_shower_energies->Add(g_shower_10);
    mg_shower_energies->Add(g_shower_20);
    mg_shower_energies->Add(g_shower_40);
    mg_shower_energies->Add(g_shower_80);
    mg_shower_energies->Add(g_shower_150);

    auto c_shower_energies = new TCanvas("c_shower_energies", "c_shower_energies", 800, 800);
    mg_shower_energies->Draw("AP");
    mg_shower_energies->SetTitle("Shower profile");

    if (transformed) mg_shower_energies->GetXaxis()->SetTitle("Layer");
    else mg_shower_energies->GetXaxis()->SetTitle("Layer");
    mg_shower_energies->GetYaxis()->SetTitle("Hits");
    mg_shower_energies->GetYaxis()->SetRangeUser(0,150);

    TLegend *leg_shower_energies;
    if(transformed == true) leg_shower_energies= new TLegend(0.45,0.65,0.65,0.85);
    else leg_shower_energies= new TLegend(0.45,0.65,0.65,0.85);
    leg_shower_energies->SetTextSize(0.035);
    leg_shower_energies->SetTextFont(42);
    leg_shower_energies->AddEntry(g_shower_2,"2 GeV","ap");
    leg_shower_energies->AddEntry(g_shower_4,"4 GeV","ap");
    leg_shower_energies->AddEntry(g_shower_6,"6 GeV","ap");
    leg_shower_energies->AddEntry(g_shower_10,"10 GeV","ap");
    leg_shower_energies->AddEntry(g_shower_20,"20 GeV","ap");
    leg_shower_energies->AddEntry(g_shower_40,"40 GeV","ap");
    leg_shower_energies->AddEntry(g_shower_80,"80 GeV","ap");
    leg_shower_energies->AddEntry(g_shower_150,"150 GeV","ap");
    leg_shower_energies->SetFillColor(0);
    leg_shower_energies->SetLineColor(0);
    leg_shower_energies->SetShadowColor(0);
    leg_shower_energies->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_energies->SaveAs("hits_shower_energies_"+savetrans+".eps");
    }


    // Shower profile normalized
    // Example layer 3                                                                                                                                          
    TGraphErrors *g_nhits_layer_w_3 = new TGraphErrors(*energyaxis, *mu_nhit_layer_w_3, *zeros, *zeros);
    graph_setup_add(g_nhits_layer_w_3, "Hits in layer 3 (weighted)", kBlue);

    TMultiGraph *mg_nhits_layer_w_3 = new TMultiGraph();
    mg_nhits_layer_w_3->Add(g_nhits_layer_w_3);

    auto c_layer_w_3 = new TCanvas("c_layer_w_3", "c_layer_w_3", 800, 800);
    mg_nhits_layer_w_3->Draw("AP");
    mg_nhits_layer_w_3->SetTitle("Hits in layer 3 (weighted)");

    if (transformed) mg_nhits_layer_w_3->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_nhits_layer_w_3->GetXaxis()->SetTitle("E[GeV]");
    mg_nhits_layer_w_3->GetYaxis()->SetTitle("Hits");
    mg_nhits_layer_w_3->GetYaxis()->SetRangeUser(0.,0.5);

    TLegend *leg_layer_w_3;
    if(transformed == true) leg_layer_w_3= new TLegend(0.45,0.65,0.65,0.85);
    else leg_layer_w_3= new TLegend(0.45,0.65,0.65,0.85);
    leg_layer_w_3->SetTextSize(0.035);
    leg_layer_w_3->SetTextFont(42);
    leg_layer_w_3->AddEntry(g_mol,"Hits in layer 3 (weighted)","ap");
    leg_layer_w_3->SetFillColor(0);
    leg_layer_w_3->SetLineColor(0);
    leg_layer_w_3->SetShadowColor(0);
    leg_layer_w_3->Draw();
    addCaliceLogo();
    if(save == true){
      c_layer_3->SaveAs("hits_layer_w_3_"+savetrans+".eps");
    }

    // Shower different energies weighted                                                                                                                                                            
    TGraphErrors *g_shower_w_2 = new TGraphErrors(*layeraxis, *shower_profile_w_2 , *zeros, *zeros);
    graph_setup_add(g_shower_w_2, "Shower profile w (2 GeV)", kBlack);
    TGraphErrors *g_shower_w_4 = new TGraphErrors(*layeraxis, *shower_profile_w_4 , *zeros, *zeros);
    graph_setup_add(g_shower_w_4, "Shower profile w (4 GeV)", kBlue+1);
    TGraphErrors *g_shower_w_6 = new TGraphErrors(*layeraxis, *shower_profile_w_6 , *zeros, *zeros);
    graph_setup_add(g_shower_w_6, "Shower profile w (6 GeV)", kBlue);
    TGraphErrors *g_shower_w_10 = new TGraphErrors(*layeraxis, *shower_profile_w_10 , *zeros, *zeros);
    graph_setup_add(g_shower_w_10, "Shower profile w (10 GeV)", kAzure-4);
    TGraphErrors *g_shower_w_20 = new TGraphErrors(*layeraxis, *shower_profile_w_20 , *zeros, *zeros);
    graph_setup_add(g_shower_w_20, "Shower profile w (20 GeV)", kAzure+8);
    TGraphErrors *g_shower_w_40 = new TGraphErrors(*layeraxis, *shower_profile_w_40 , *zeros, *zeros);
    graph_setup_add(g_shower_w_40, "Shower profile w (40 GeV)", kRed-7);
    TGraphErrors *g_shower_w_80 = new TGraphErrors(*layeraxis, *shower_profile_w_80 , *zeros, *zeros);
    graph_setup_add(g_shower_w_80, "Shower profile w (80 GeV)", kRed);
    TGraphErrors *g_shower_w_150 = new TGraphErrors(*layeraxis, *shower_profile_w_150 , *zeros, *zeros);
    graph_setup_add(g_shower_w_150, "Shower profile w (150 GeV)", kRed-3);

    TMultiGraph *mg_shower_energies_w = new TMultiGraph();
    mg_shower_energies_w->Add(g_shower_w_2);
    mg_shower_energies_w->Add(g_shower_w_4);
    mg_shower_energies_w->Add(g_shower_w_6);
    mg_shower_energies_w->Add(g_shower_w_10);
    mg_shower_energies_w->Add(g_shower_w_20);
    mg_shower_energies_w->Add(g_shower_w_40);
    mg_shower_energies_w->Add(g_shower_w_80);
    mg_shower_energies_w->Add(g_shower_w_150);

    auto c_shower_energies_w = new TCanvas("c_shower_energies_w", "c_shower_energies_w", 800, 800);
    mg_shower_energies_w->Draw("AP");
    mg_shower_energies_w->SetTitle("Shower profile (normalized)");

    if (transformed) mg_shower_energies_w->GetXaxis()->SetTitle("Layer");
    else mg_shower_energies_w->GetXaxis()->SetTitle("Layer");
    mg_shower_energies_w->GetYaxis()->SetTitle("Hits (normalized)");
    mg_shower_energies_w->GetYaxis()->SetRangeUser(0.,0.5);

    TLegend *leg_shower_energies_w;
    if(transformed == true) leg_shower_energies_w= new TLegend(0.45,0.65,0.65,0.85);
    else leg_shower_energies_w= new TLegend(0.45,0.65,0.65,0.85);
    leg_shower_energies_w->SetTextSize(0.035);
    leg_shower_energies_w->SetTextFont(42);
    leg_shower_energies_w->AddEntry(g_shower_w_2,"2 GeV","ap");
    leg_shower_energies_w->AddEntry(g_shower_w_4,"4 GeV","ap");
    leg_shower_energies_w->AddEntry(g_shower_w_6,"6 GeV","ap");
    leg_shower_energies_w->AddEntry(g_shower_w_10,"10 GeV","ap");
    leg_shower_energies_w->AddEntry(g_shower_w_20,"20 GeV","ap");
    leg_shower_energies_w->AddEntry(g_shower_w_40,"40 GeV","ap");
    leg_shower_energies_w->AddEntry(g_shower_w_80,"80 GeV","ap");
    leg_shower_energies_w->AddEntry(g_shower_w_150,"150 GeV","ap");
    leg_shower_energies_w->SetFillColor(0);
    leg_shower_energies_w->SetLineColor(0);
    leg_shower_energies_w->SetShadowColor(0);
    leg_shower_energies_w->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_energies->SaveAs("hits_shower_energies_w_"+savetrans+".eps");
    }

    // Shower max, start and end
    // Shower nhit max
    TGraphErrors *g_shower_nhit_max = new TGraphErrors(*energyaxis, *shower_nhit_max, *zeros, *zeros);
    graph_setup_add(g_shower_nhit_max, "Shower maximum (n hits)", kBlue);

    TMultiGraph *mg_shower_nhit_max = new TMultiGraph();
    mg_shower_nhit_max->Add(g_shower_nhit_max);

    auto c_shower_nhit_max = new TCanvas("c_shower_nhit_max", "c_shower_nhit_max", 800, 800);
    mg_shower_nhit_max->Draw("AP");
    mg_shower_nhit_max->SetTitle("Shower max (n hits)");

    if (transformed) mg_shower_nhit_max->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_nhit_max->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_nhit_max->GetYaxis()->SetTitle("Layer");
    mg_shower_nhit_max->GetYaxis()->SetRangeUser(0,50);

    TLegend *leg_shower_nhit_max;
    if(transformed == true) leg_shower_nhit_max= new TLegend(0.45,0.65,0.65,0.85);
    else leg_shower_nhit_max= new TLegend(0.45,0.65,0.65,0.85);
    leg_shower_nhit_max->SetTextSize(0.035);
    leg_shower_nhit_max->SetTextFont(42);
    leg_shower_nhit_max->AddEntry(g_shower_nhit_max,"Shower max. (n hits)","ap");
    leg_shower_nhit_max->SetFillColor(0);
    leg_shower_nhit_max->SetLineColor(0);
    leg_shower_nhit_max->SetShadowColor(0);
    leg_shower_nhit_max->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_nhit_max->SaveAs("hits_shower_nhit_max_"+savetrans+".eps");
    }

    // Shower nhit average
    TGraphErrors *g_shower_nhit_average = new TGraphErrors(*energyaxis, *shower_nhit_average, *zeros, *zeros);
    graph_setup_add(g_shower_nhit_average, "Shower average (n hits)", kBlue);

    TMultiGraph *mg_shower_nhit_average = new TMultiGraph();
    mg_shower_nhit_average->Add(g_shower_nhit_average);

    auto c_shower_nhit_average = new TCanvas("c_shower_nhit_average", "c_shower_nhit_average", 800, 800);
    mg_shower_nhit_average->Draw("AP");
    mg_shower_nhit_average->SetTitle("Shower average (n hits)");

    if (transformed) mg_shower_nhit_average->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_nhit_average->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_nhit_average->GetYaxis()->SetTitle("Layer");
    mg_shower_nhit_average->GetYaxis()->SetRangeUser(0,50);

    TLegend *leg_shower_nhit_average;
    if(transformed == true) leg_shower_nhit_average= new TLegend(0.45,0.65,0.65,0.85);
    else leg_shower_nhit_average= new TLegend(0.45,0.65,0.65,0.85);
    leg_shower_nhit_average->SetTextSize(0.035);
    leg_shower_nhit_average->SetTextFont(42);
    leg_shower_nhit_average->AddEntry(g_shower_nhit_average,"Shower average (n hits)","ap");
    leg_shower_nhit_average->SetFillColor(0);
    leg_shower_nhit_average->SetLineColor(0);
    leg_shower_nhit_average->SetShadowColor(0);
    leg_shower_nhit_average->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_nhit_average->SaveAs("hits_shower_nhit_average_"+savetrans+".eps");
    }

    // Shower nhit max layer
    TGraphErrors *g_shower_nhit_max_layer = new TGraphErrors(*energyaxis, *shower_nhit_max_layer, *zeros, *zeros);
    graph_setup_add(g_shower_nhit_max_layer, "Shower maximum (n hits)", kBlue);

    TMultiGraph *mg_shower_nhit_max_layer = new TMultiGraph();
    mg_shower_nhit_max_layer->Add(g_shower_nhit_max_layer);

    auto c_shower_nhit_max_layer = new TCanvas("c_shower_nhit_max_layer", "c_shower_nhit_max_layer", 800, 800);
    mg_shower_nhit_max_layer->Draw("AP");
    mg_shower_nhit_max_layer->SetTitle("N hits shower max layer (avg.)");

    if (transformed) mg_shower_nhit_max_layer->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_nhit_max_layer->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_nhit_max_layer->GetYaxis()->SetTitle("Layer");
    mg_shower_nhit_max_layer->GetYaxis()->SetRangeUser(-1,15);

    TLegend *leg_shower_nhit_max_layer;
    if(transformed == true) leg_shower_nhit_max_layer= new TLegend(0.45,0.65,0.65,0.85);
    else leg_shower_nhit_max_layer= new TLegend(0.45,0.65,0.65,0.85);
    leg_shower_nhit_max_layer->SetTextSize(0.035);
    leg_shower_nhit_max_layer->SetTextFont(42);
    leg_shower_nhit_max_layer->AddEntry(g_shower_nhit_max_layer,"Shower max layer (n hits)","ap");
    leg_shower_nhit_max_layer->SetFillColor(0);
    leg_shower_nhit_max_layer->SetLineColor(0);
    leg_shower_nhit_max_layer->SetShadowColor(0);
    leg_shower_nhit_max_layer->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_nhit_max_layer->SaveAs("hits_shower_nhit_max_layer_"+savetrans+".eps");
    }

    // Shower nhit start layer
    TGraphErrors *g_shower_nhit_start_layer = new TGraphErrors(*energyaxis, *shower_nhit_start_layer, *zeros, *zeros);
    graph_setup_add(g_shower_nhit_start_layer, "Shower startimum (n hits)", kBlue);

    TMultiGraph *mg_shower_nhit_start_layer = new TMultiGraph();
    mg_shower_nhit_start_layer->Add(g_shower_nhit_start_layer);

    auto c_shower_nhit_start_layer = new TCanvas("c_shower_nhit_start_layer", "c_shower_nhit_start_layer", 800, 800);
    mg_shower_nhit_start_layer->Draw("AP");
    mg_shower_nhit_start_layer->SetTitle("N hits shower start layer (avg.)");

    if (transformed) mg_shower_nhit_start_layer->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_nhit_start_layer->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_nhit_start_layer->GetYaxis()->SetTitle("Layer");
    mg_shower_nhit_start_layer->GetYaxis()->SetRangeUser(-1,15);

    TLegend *leg_shower_nhit_start_layer;
    if(transformed == true) leg_shower_nhit_start_layer= new TLegend(0.45,0.65,0.65,0.85);
    else leg_shower_nhit_start_layer= new TLegend(0.45,0.65,0.65,0.85);
    leg_shower_nhit_start_layer->SetTextSize(0.035);
    leg_shower_nhit_start_layer->SetTextFont(42);
    leg_shower_nhit_start_layer->AddEntry(g_shower_nhit_start_layer,"Shower start layer (n hits)","ap");
    leg_shower_nhit_start_layer->SetFillColor(0);
    leg_shower_nhit_start_layer->SetLineColor(0);
    leg_shower_nhit_start_layer->SetShadowColor(0);
    leg_shower_nhit_start_layer->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_nhit_start_layer->SaveAs("hits_shower_nhit_start_layer_"+savetrans+".eps");
    }

    // Shower nhit end layer                                                                                                                                                                                                            
    TGraphErrors *g_shower_nhit_end_layer = new TGraphErrors(*energyaxis, *shower_nhit_end_layer, *zeros, *zeros);
    graph_setup_add(g_shower_nhit_end_layer, "Shower endimum (n hits)", kBlue);

    TMultiGraph *mg_shower_nhit_end_layer = new TMultiGraph();
    mg_shower_nhit_end_layer->Add(g_shower_nhit_end_layer);

    auto c_shower_nhit_end_layer = new TCanvas("c_shower_nhit_end_layer", "c_shower_nhit_end_layer", 800, 800);
    mg_shower_nhit_end_layer->Draw("AP");
    mg_shower_nhit_end_layer->SetTitle("N hits shower end layer (avg.)");

    if (transformed) mg_shower_nhit_end_layer->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_nhit_end_layer->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_nhit_end_layer->GetYaxis()->SetTitle("Layer");
    mg_shower_nhit_end_layer->GetYaxis()->SetRangeUser(-1,15);

    TLegend *leg_shower_nhit_end_layer;
    if(transformed == true) leg_shower_nhit_end_layer= new TLegend(0.45,0.65,0.65,0.85);
    else leg_shower_nhit_end_layer= new TLegend(0.45,0.65,0.65,0.85);
    leg_shower_nhit_end_layer->SetTextSize(0.035);
    leg_shower_nhit_end_layer->SetTextFont(42);
    leg_shower_nhit_end_layer->AddEntry(g_shower_nhit_end_layer,"Shower end layer (n hits)","ap");
    leg_shower_nhit_end_layer->SetFillColor(0);
    leg_shower_nhit_end_layer->SetLineColor(0);
    leg_shower_nhit_end_layer->SetShadowColor(0);
    leg_shower_nhit_end_layer->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_nhit_end_layer->SaveAs("hits_shower_nhit_end_layer_"+savetrans+".eps");
    }


}
