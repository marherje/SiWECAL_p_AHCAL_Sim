#define N_ENERGIES 18
#define N_ECAL_LAYERS 15
// We add the -1 layer

void addCaliceLogo(bool WIP = true){
  TImage *img = TImage::Open("../style/CALICELogo_18pc.png");
  img->SetConstRatio(kTRUE);
  img->SetImageCompression(0);
  TPad *p1 = new TPad("img", "img", 0.835, 0.92, 1.0, 1.0);
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

void res_mol_plots_muerrors(string particle, bool transformed = true, bool save = false){
    
    TString filename = "results_folder/resolution_"+particle+"_result.root";
    TFile *file = new TFile(filename, "read");
    // TVectorD * energies = file->GetObject("energies");
    TString savetrans;
    if(transformed == false) savetrans="Linear_"+particle;
    else savetrans="InvSq_"+particle;

    TString partstring(particle);
    
    TVectorD *energies = new TVectorD(N_ENERGIES);
    TVectorD *energies_tr  = new TVectorD(N_ENERGIES);
    
    TVectorD *res_nhit  = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit = new TVectorD(N_ENERGIES);
    TVectorD *res_sume = new TVectorD(N_ENERGIES);
    TVectorD *mu_sume = new TVectorD(N_ENERGIES);
    TVectorD *sig_sume = new TVectorD(N_ENERGIES);
    TVectorD *res_weight = new TVectorD(N_ENERGIES);
    TVectorD *mu_weight = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight = new TVectorD(N_ENERGIES);
    
    TVectorD *mu_mol = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_mol = new TVectorD(N_ENERGIES);
    TVectorD *sig_mol = new TVectorD(N_ENERGIES);

    TVectorD *mu_barycenter_x = new TVectorD(N_ENERGIES);
    TVectorD *mu_barycenter_y = new TVectorD(N_ENERGIES);
    TVectorD *mu_barycenter_z = new TVectorD(N_ENERGIES);
    TVectorD *sig_barycenter_x = new TVectorD(N_ENERGIES);
    TVectorD *sig_barycenter_y = new TVectorD(N_ENERGIES);
    TVectorD *sig_barycenter_z = new TVectorD(N_ENERGIES);

    TVectorD *mu_nhit_layer_0 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_nhit_layer_1 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_nhit_layer_2 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_nhit_layer_3 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_4 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_nhit_layer_5 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_nhit_layer_6 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_7 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_8 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_9 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_nhit_layer_10 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_nhit_layer_11 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_nhit_layer_12 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_nhit_layer_13 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_nhit_layer_14 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_error_nhit_layer_0 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_1 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_2 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_3 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_4 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_5 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_6 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_7 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_8 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_9 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_10 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_11 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_12 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_13 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_14 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_0 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_1 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_2 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_3 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_4 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_5 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_6 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_7 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_8 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_9 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_10 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_11 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_12 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_13 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_14 = new TVectorD(N_ENERGIES);
    
    TVectorD *mu_nhit_layer_n_0 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_1 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_2 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_3 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_4 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_5 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_6 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_7 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_8 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_9 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_10 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_11 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_12 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_13 = new TVectorD(N_ENERGIES);
    TVectorD *mu_nhit_layer_n_14 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_0 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_1 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_2 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_3 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_4 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_5 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_6 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_7 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_8 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_9 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_10 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_11 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_12 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_13 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_nhit_layer_n_14 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_0 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_1 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_2 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_3 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_4 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_5 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_6 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_7 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_8 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_9 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_10 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_11 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_12 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_13 = new TVectorD(N_ENERGIES);
    TVectorD *sig_nhit_layer_n_14 = new TVectorD(N_ENERGIES);
    TVectorD *shower_nhit_max_layer = new TVectorD(N_ENERGIES); 
    TVectorD *shower_nhit_start_layer = new TVectorD(N_ENERGIES); 
    TVectorD *shower_nhit_end_layer = new TVectorD(N_ENERGIES); 
    TVectorD *shower_nhit_max = new TVectorD(N_ENERGIES); 
    TVectorD *shower_nhit_average = new TVectorD(N_ENERGIES);

    TVectorD *mu_weight_layer_0 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_1 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_2 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_3 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_4 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_5 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_6 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_7 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_8 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_9 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_10 = new TVectorD(N_ENERGIES);
    TVectorD *mu_weight_layer_11 = new TVectorD(N_ENERGIES);
    TVectorD *mu_weight_layer_12 = new TVectorD(N_ENERGIES);
    TVectorD *mu_weight_layer_13 = new TVectorD(N_ENERGIES);
    TVectorD *mu_weight_layer_14 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_0 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_1 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_2 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_3 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_4 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_5 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_6 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_7 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_8 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_9 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_10 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_11 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_12 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_13 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_14 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_0 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_1 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_2 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_3 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_4 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_5 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_6 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_7 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_8 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_9 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_10 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_11 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_12 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_13 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_14 = new TVectorD(N_ENERGIES);

    TVectorD *mu_weight_layer_n_0 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_n_1 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_n_2 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_n_3 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_n_4 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_n_5 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_n_6 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_n_7 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_n_8 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_n_9 = new TVectorD(N_ENERGIES); 
    TVectorD *mu_weight_layer_n_10 = new TVectorD(N_ENERGIES);
    TVectorD *mu_weight_layer_n_11 = new TVectorD(N_ENERGIES);
    TVectorD *mu_weight_layer_n_12 = new TVectorD(N_ENERGIES);
    TVectorD *mu_weight_layer_n_13 = new TVectorD(N_ENERGIES);
    TVectorD *mu_weight_layer_n_14 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_0 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_1 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_2 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_3 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_4 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_5 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_6 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_7 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_8 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_9 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_10 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_11 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_12 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_13 = new TVectorD(N_ENERGIES);
    TVectorD *mu_error_weight_layer_n_14 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_0 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_1 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_2 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_3 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_4 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_5 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_6 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_7 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_8 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_9 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_10 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_11 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_12 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_13 = new TVectorD(N_ENERGIES);
    TVectorD *sig_weight_layer_n_14 = new TVectorD(N_ENERGIES);
    
    TVectorD *shower_weight_max_layer = new TVectorD(N_ENERGIES); 
    TVectorD *shower_weight_start_layer = new TVectorD(N_ENERGIES); 
    TVectorD *shower_weight_end_layer = new TVectorD(N_ENERGIES); 
    TVectorD *shower_weight_max  = new TVectorD(N_ENERGIES); 
    TVectorD *shower_weight_average = new TVectorD(N_ENERGIES);


    // New variables
    // Index    0 1 2 3  4  5  6  7  8  9 10 11 12  13  14  15  16  17
    // Energies 2 4 6 8 10 20 30 40 50 60 70 80 90 100 125 150 175 200
    // Shower nhit
    double showernhitvalues[N_ENERGIES][N_ECAL_LAYERS];
    double showernhitvalues_n[N_ENERGIES][N_ECAL_LAYERS];

    double showernhiterrors[N_ENERGIES][N_ECAL_LAYERS];
    double showernhiterrors_n[N_ENERGIES][N_ECAL_LAYERS];

    TVectorD *layeraxis = new TVectorD(N_ECAL_LAYERS);
    double layeraxisvalues[N_ECAL_LAYERS] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.};
    TVectorD layeraxisref(N_ECAL_LAYERS, layeraxisvalues);
    layeraxis = &layeraxisref;
    
    TVectorD *shower_nhit_profile_2 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_2[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_4 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_4[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_6 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_6[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_10 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_10[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_20 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_20[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_40 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_40[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_80 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_80[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_150 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_150[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_2 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_2[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_4 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_4[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_6 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_6[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_10 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_10[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_20 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_20[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_40 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_40[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_80 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_80[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_150 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_150[N_ECAL_LAYERS];
    
    TVectorD *shower_nhit_profile_n_2 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_n_2[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_n_4 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_n_4[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_n_6 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_n_6[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_n_10 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_n_10[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_n_20 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_n_20[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_n_40 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_n_40[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_n_80 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_n_80[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_n_150 = new TVectorD(N_ECAL_LAYERS);
    double showernhitvalues_n_150[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_n_2 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_n_2[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_n_4 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_n_4[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_n_6 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_n_6[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_n_10 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_n_10[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_n_20 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_n_20[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_n_40 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_n_40[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_n_80 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_n_80[N_ECAL_LAYERS];
    TVectorD *shower_nhit_profile_errors_n_150 = new TVectorD(N_ECAL_LAYERS);
    double showernhiterrors_n_150[N_ECAL_LAYERS];

    TVectorD *zeros = new TVectorD(N_ENERGIES);
    double zerovalues[N_ENERGIES];
    for (int i=0; i<N_ENERGIES; i++) zerovalues[i] = 0.;
    TVectorD zerovectorref(N_ENERGIES, zerovalues);
    zeros = &zerovectorref;

    // +1 fr the -1 layer
    TVectorD *zeroslayer = new TVectorD(N_ECAL_LAYERS);
    double zerovalueslayer[N_ECAL_LAYERS];
    for(int i=0; i<N_ECAL_LAYERS;i++) zerovalueslayer[i] = 0.;
    TVectorD zerovectorreflayer(N_ECAL_LAYERS, zerovalueslayer);
    zeroslayer = &zerovectorreflayer;
    
    // Shower weighted energy
    double showerweightvalues[N_ENERGIES][N_ECAL_LAYERS];
    double showerweightvalues_n[N_ENERGIES][N_ECAL_LAYERS];

    double showerweighterrors[N_ENERGIES][N_ECAL_LAYERS];
    double showerweighterrors_n[N_ENERGIES][N_ECAL_LAYERS];

    TVectorD *shower_weight_profile_2 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_2[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_4 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_4[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_6 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_6[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_10 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_10[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_20 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_20[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_40 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_40[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_80 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_80[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_150 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_150[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_2 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_2[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_4 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_4[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_6 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_6[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_10 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_10[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_20 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_20[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_40 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_40[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_80 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_80[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_150 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_150[N_ECAL_LAYERS];

    TVectorD *shower_weight_profile_n_2 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_n_2[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_n_4 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_n_4[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_n_6 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_n_6[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_n_10 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_n_10[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_n_20 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_n_20[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_n_40 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_n_40[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_n_80 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_n_80[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_n_150 = new TVectorD(N_ECAL_LAYERS);
    double showerweightvalues_n_150[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_n_2 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_n_2[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_n_4 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_n_4[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_n_6 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_n_6[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_n_10 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_n_10[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_n_20 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_n_20[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_n_40 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_n_40[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_n_80 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_n_80[N_ECAL_LAYERS];
    TVectorD *shower_weight_profile_errors_n_150 = new TVectorD(N_ECAL_LAYERS);
    double showerweighterrors_n_150[N_ECAL_LAYERS];

    // Load objects
    
    for(int i=1; i<N_ENERGIES+1; i++) {
    string icycle = to_string(i);
    cout<<"cycle: "+icycle<<endl;
    TVectorD *energies_cycle = (TVectorD*)(file->FindObjectAny(("energies;"+icycle).c_str()));
    TVectorD *energies_tr_cycle = (TVectorD*)(file->FindObjectAny(("energies_tr;"+icycle).c_str()));
    TVectorD *res_nhit_cycle = (TVectorD*)(file->FindObjectAny(("res_nhit;"+icycle).c_str()));
    TVectorD *res_sume_cycle = (TVectorD*)(file->FindObjectAny(("res_sume;"+icycle).c_str()));
    TVectorD *res_weight_cycle = (TVectorD*)(file->FindObjectAny(("res_weight;"+icycle).c_str()));

    TVectorD *mu_nhit_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit;"+icycle).c_str()));
    TVectorD *mu_sume_cycle = (TVectorD*)(file->FindObjectAny(("mu_sume;"+icycle).c_str()));
    TVectorD *mu_weight_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight;"+icycle).c_str()));

    TVectorD *sig_nhit_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit;"+icycle).c_str()));
    TVectorD *sig_sume_cycle = (TVectorD*)(file->FindObjectAny(("sig_sume;"+icycle).c_str()));
    TVectorD *sig_weight_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight;"+icycle).c_str()));

    TVectorD *mu_mol_cycle = (TVectorD*)(file->FindObjectAny(("mu_mol;"+icycle).c_str()));
    TVectorD *sig_mol_cycle = (TVectorD*)(file->FindObjectAny(("sig_mol;"+icycle).c_str()));

    TVectorD *mu_barycenter_x_cycle = (TVectorD*)(file->FindObjectAny(("mu_barycenter_x;"+icycle).c_str()));
    TVectorD *mu_barycenter_y_cycle = (TVectorD*)(file->FindObjectAny(("mu_barycenter_y;"+icycle).c_str()));
    TVectorD *mu_barycenter_z_cycle = (TVectorD*)(file->FindObjectAny(("mu_barycenter_z;"+icycle).c_str()));
    TVectorD *sig_barycenter_x_cycle = (TVectorD*)(file->FindObjectAny(("sig_barycenter_x;"+icycle).c_str()));
    TVectorD *sig_barycenter_y_cycle = (TVectorD*)(file->FindObjectAny(("sig_barycenter_y;"+icycle).c_str()));
    TVectorD *sig_barycenter_z_cycle = (TVectorD*)(file->FindObjectAny(("sig_barycenter_z;"+icycle).c_str()));

    TVectorD *mu_nhit_layer_0_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_0;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_1_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_1;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_2_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_2;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_3_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_3;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_4_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_4;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_5_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_5;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_6_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_6;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_7_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_7;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_8_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_8;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_9_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_9;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_10_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_10;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_11_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_11;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_12_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_12;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_13_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_13;"+icycle).c_str()));
    TVectorD *mu_nhit_layer_14_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_14;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_0_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_0;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_1_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_1;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_2_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_2;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_3_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_3;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_4_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_4;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_5_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_5;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_6_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_6;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_7_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_7;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_8_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_8;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_9_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_9;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_10_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_10;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_11_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_11;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_12_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_12;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_13_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_13;"+icycle).c_str()));
    TVectorD *mu_error_nhit_layer_14_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_14;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_0_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_0;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_1_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_1;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_2_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_2;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_3_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_3;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_4_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_4;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_5_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_5;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_6_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_6;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_7_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_7;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_8_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_8;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_9_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_9;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_10_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_10;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_11_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_11;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_12_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_12;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_13_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_13;"+icycle).c_str()));
    TVectorD *sig_nhit_layer_14_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_14;"+icycle).c_str()));

    TVectorD * mu_nhit_layer_n_0_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_0;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_1_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_1;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_2_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_2;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_3_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_3;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_4_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_4;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_5_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_5;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_6_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_6;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_7_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_7;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_8_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_8;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_9_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_9;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_10_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_10;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_11_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_11;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_12_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_12;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_13_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_13;"+icycle).c_str()));
    TVectorD * mu_nhit_layer_n_14_cycle = (TVectorD*)(file->FindObjectAny(("mu_nhit_layer_n_14;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_0_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_0;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_1_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_1;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_2_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_2;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_3_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_3;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_4_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_4;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_5_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_5;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_6_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_6;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_7_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_7;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_8_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_8;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_9_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_9;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_10_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_10;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_11_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_11;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_12_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_12;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_13_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_13;"+icycle).c_str()));
    TVectorD * mu_error_nhit_layer_n_14_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_nhit_layer_n_14;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_0_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_0;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_1_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_1;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_2_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_2;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_3_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_3;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_4_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_4;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_5_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_5;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_6_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_6;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_7_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_7;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_8_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_8;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_9_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_9;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_10_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_10;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_11_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_11;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_12_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_12;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_13_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_13;"+icycle).c_str()));
    TVectorD * sig_nhit_layer_n_14_cycle = (TVectorD*)(file->FindObjectAny(("sig_nhit_layer_n_14;"+icycle).c_str()));
    TVectorD * shower_nhit_max_layer_cycle = (TVectorD*)(file->FindObjectAny(("shower_nhit_max_layer;"+icycle).c_str()));
    TVectorD * shower_nhit_start_layer_cycle = (TVectorD*)(file->FindObjectAny(("shower_nhit_start_layer;"+icycle).c_str()));
    TVectorD * shower_nhit_end_layer_cycle = (TVectorD*)(file->FindObjectAny(("shower_nhit_end_layer;"+icycle).c_str()));
    TVectorD * shower_nhit_max_cycle = (TVectorD*)(file->FindObjectAny(("shower_nhit_max;"+icycle).c_str()));
    TVectorD * shower_nhit_average_cycle = (TVectorD*)(file->FindObjectAny(("shower_nhit_average;"+icycle).c_str()));

    //weighted energy
    TVectorD * mu_weight_layer_0_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_0;"+icycle).c_str()));
    TVectorD * mu_weight_layer_1_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_1;"+icycle).c_str()));
    TVectorD * mu_weight_layer_2_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_2;"+icycle).c_str()));
    TVectorD * mu_weight_layer_3_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_3;"+icycle).c_str()));
    TVectorD * mu_weight_layer_4_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_4;"+icycle).c_str()));
    TVectorD * mu_weight_layer_5_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_5;"+icycle).c_str()));
    TVectorD * mu_weight_layer_6_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_6;"+icycle).c_str()));
    TVectorD * mu_weight_layer_7_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_7;"+icycle).c_str()));
    TVectorD * mu_weight_layer_8_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_8;"+icycle).c_str()));
    TVectorD * mu_weight_layer_9_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_9;"+icycle).c_str()));
    TVectorD * mu_weight_layer_10_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_10;"+icycle).c_str()));
    TVectorD * mu_weight_layer_11_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_11;"+icycle).c_str()));
    TVectorD * mu_weight_layer_12_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_12;"+icycle).c_str()));
    TVectorD * mu_weight_layer_13_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_13;"+icycle).c_str()));
    TVectorD * mu_weight_layer_14_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_14;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_0_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_0;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_1_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_1;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_2_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_2;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_3_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_3;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_4_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_4;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_5_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_5;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_6_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_6;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_7_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_7;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_8_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_8;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_9_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_9;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_10_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_10;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_11_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_11;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_12_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_12;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_13_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_13;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_14_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_14;"+icycle).c_str()));
    TVectorD * sig_weight_layer_0_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_0;"+icycle).c_str()));
    TVectorD * sig_weight_layer_1_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_1;"+icycle).c_str()));
    TVectorD * sig_weight_layer_2_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_2;"+icycle).c_str()));
    TVectorD * sig_weight_layer_3_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_3;"+icycle).c_str()));
    TVectorD * sig_weight_layer_4_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_4;"+icycle).c_str()));
    TVectorD * sig_weight_layer_5_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_5;"+icycle).c_str()));
    TVectorD * sig_weight_layer_6_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_6;"+icycle).c_str()));
    TVectorD * sig_weight_layer_7_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_7;"+icycle).c_str()));
    TVectorD * sig_weight_layer_8_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_8;"+icycle).c_str()));
    TVectorD * sig_weight_layer_9_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_9;"+icycle).c_str()));
    TVectorD * sig_weight_layer_10_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_10;"+icycle).c_str()));
    TVectorD * sig_weight_layer_11_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_11;"+icycle).c_str()));
    TVectorD * sig_weight_layer_12_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_12;"+icycle).c_str()));
    TVectorD * sig_weight_layer_13_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_13;"+icycle).c_str()));
    TVectorD * sig_weight_layer_14_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_14;"+icycle).c_str()));

    TVectorD * mu_weight_layer_n_0_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_0;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_1_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_1;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_2_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_2;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_3_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_3;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_4_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_4;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_5_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_5;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_6_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_6;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_7_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_7;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_8_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_8;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_9_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_9;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_10_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_10;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_11_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_11;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_12_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_12;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_13_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_13;"+icycle).c_str()));
    TVectorD * mu_weight_layer_n_14_cycle = (TVectorD*)(file->FindObjectAny(("mu_weight_layer_n_14;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_0_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_0;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_1_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_1;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_2_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_2;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_3_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_3;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_4_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_4;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_5_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_5;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_6_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_6;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_7_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_7;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_8_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_8;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_9_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_9;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_10_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_10;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_11_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_11;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_12_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_12;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_13_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_13;"+icycle).c_str()));
    TVectorD * mu_error_weight_layer_n_14_cycle = (TVectorD*)(file->FindObjectAny(("mu_error_weight_layer_n_14;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_0_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_0;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_1_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_1;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_2_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_2;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_3_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_3;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_4_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_4;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_5_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_5;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_6_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_6;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_7_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_7;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_8_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_8;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_9_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_9;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_10_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_10;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_11_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_11;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_12_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_12;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_13_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_13;"+icycle).c_str()));
    TVectorD * sig_weight_layer_n_14_cycle = (TVectorD*)(file->FindObjectAny(("sig_weight_layer_n_14;"+icycle).c_str()));

    TVectorD * shower_weight_max_layer_cycle = (TVectorD*)(file->FindObjectAny(("shower_weight_max_layer;"+icycle).c_str()));
    TVectorD * shower_weight_start_layer_cycle = (TVectorD*)(file->FindObjectAny(("shower_weight_start_layer;"+icycle).c_str()));
    TVectorD * shower_weight_end_layer_cycle = (TVectorD*)(file->FindObjectAny(("shower_weight_end_layer;"+icycle).c_str()));
    TVectorD * shower_weight_max_cycle = (TVectorD*)(file->FindObjectAny(("shower_weight_max;"+icycle).c_str()));
    TVectorD * shower_weight_average_cycle = (TVectorD*)(file->FindObjectAny(("shower_weight_average;"+icycle).c_str()));

    ((*energies))[i-1] = ((*energies_cycle))[0];
    ((*energies_tr))[i-1] = ((*energies_tr_cycle))[0];
    ((*res_nhit))[i-1] = ((*res_nhit_cycle))[0];
    ((*res_sume))[i-1] = ((*res_sume_cycle))[0];
    ((*res_weight))[i-1] = ((*res_weight_cycle))[0];

    ((*mu_nhit))[i-1] = ((*mu_nhit_cycle))[0];
    ((*mu_sume))[i-1] = ((*mu_sume_cycle))[0];
    ((*mu_weight))[i-1] = ((*mu_weight_cycle))[0];

    ((*sig_nhit))[i-1] = ((*sig_nhit_cycle))[0];
    ((*sig_sume))[i-1] = ((*sig_sume_cycle))[0];
    ((*sig_weight))[i-1] = ((*sig_weight_cycle))[0];

    ((*mu_mol))[i-1] = ((*mu_mol_cycle))[0];
    ((*sig_mol))[i-1] = ((*sig_mol_cycle))[0];

    ((*mu_barycenter_x))[i-1] = ((*mu_barycenter_x_cycle))[0];
    ((*mu_barycenter_y))[i-1] = ((*mu_barycenter_y_cycle))[0];
    ((*mu_barycenter_z))[i-1] = ((*mu_barycenter_z_cycle))[0];
    ((*sig_barycenter_x))[i-1] = ((*sig_barycenter_x_cycle))[0];
    ((*sig_barycenter_y))[i-1] = ((*sig_barycenter_y_cycle))[0];
    ((*sig_barycenter_z))[i-1] = ((*sig_barycenter_z_cycle))[0];

    ((*mu_nhit_layer_0))[i-1] = ((*mu_nhit_layer_0_cycle))[0]; 
    ((*mu_nhit_layer_1))[i-1] = ((*mu_nhit_layer_1_cycle))[0];
    ((*mu_nhit_layer_2))[i-1] = ((*mu_nhit_layer_2_cycle))[0];
    ((*mu_nhit_layer_3))[i-1] = ((*mu_nhit_layer_3_cycle))[0];
    ((*mu_nhit_layer_4))[i-1] = ((*mu_nhit_layer_4_cycle))[0];
    ((*mu_nhit_layer_5))[i-1] = ((*mu_nhit_layer_5_cycle))[0];
    ((*mu_nhit_layer_6))[i-1] = ((*mu_nhit_layer_6_cycle))[0];
    ((*mu_nhit_layer_7))[i-1] = ((*mu_nhit_layer_7_cycle))[0];
    ((*mu_nhit_layer_8))[i-1] = ((*mu_nhit_layer_8_cycle))[0];
    ((*mu_nhit_layer_9))[i-1] = ((*mu_nhit_layer_9_cycle))[0];
    ((*mu_nhit_layer_10))[i-1] = ((*mu_nhit_layer_10_cycle))[0];
    ((*mu_nhit_layer_11))[i-1] = ((*mu_nhit_layer_11_cycle))[0];
    ((*mu_nhit_layer_12))[i-1] = ((*mu_nhit_layer_12_cycle))[0];
    ((*mu_nhit_layer_13))[i-1] = ((*mu_nhit_layer_13_cycle))[0];
    ((*mu_nhit_layer_14))[i-1] = ((*mu_nhit_layer_14_cycle))[0];
    ((*mu_error_nhit_layer_0))[i-1] = ((*mu_error_nhit_layer_0_cycle))[0];
    ((*mu_error_nhit_layer_1))[i-1] = ((*mu_error_nhit_layer_1_cycle))[0];
    ((*mu_error_nhit_layer_2))[i-1] = ((*mu_error_nhit_layer_2_cycle))[0];
    ((*mu_error_nhit_layer_3))[i-1] = ((*mu_error_nhit_layer_3_cycle))[0];
    ((*mu_error_nhit_layer_4))[i-1] = ((*mu_error_nhit_layer_4_cycle))[0];
    ((*mu_error_nhit_layer_5))[i-1] = ((*mu_error_nhit_layer_5_cycle))[0];
    ((*mu_error_nhit_layer_6))[i-1] = ((*mu_error_nhit_layer_6_cycle))[0];
    ((*mu_error_nhit_layer_7))[i-1] = ((*mu_error_nhit_layer_7_cycle))[0];
    ((*mu_error_nhit_layer_8))[i-1] = ((*mu_error_nhit_layer_8_cycle))[0];
    ((*mu_error_nhit_layer_9))[i-1] = ((*mu_error_nhit_layer_9_cycle))[0];
    ((*mu_error_nhit_layer_10))[i-1] = ((*mu_error_nhit_layer_10_cycle))[0];
    ((*mu_error_nhit_layer_11))[i-1] = ((*mu_error_nhit_layer_11_cycle))[0];
    ((*mu_error_nhit_layer_12))[i-1] = ((*mu_error_nhit_layer_12_cycle))[0];
    ((*mu_error_nhit_layer_13))[i-1] = ((*mu_error_nhit_layer_13_cycle))[0];
    ((*mu_error_nhit_layer_14))[i-1] = ((*mu_error_nhit_layer_14_cycle))[0];
    ((*sig_nhit_layer_0))[i-1] = ((*sig_nhit_layer_0_cycle))[0];
    ((*sig_nhit_layer_1))[i-1] = ((*sig_nhit_layer_1_cycle))[0];
    ((*sig_nhit_layer_2))[i-1] = ((*sig_nhit_layer_2_cycle))[0];
    ((*sig_nhit_layer_3))[i-1] = ((*sig_nhit_layer_3_cycle))[0];
    ((*sig_nhit_layer_4))[i-1] = ((*sig_nhit_layer_4_cycle))[0];
    ((*sig_nhit_layer_5))[i-1] = ((*sig_nhit_layer_5_cycle))[0];
    ((*sig_nhit_layer_6))[i-1] = ((*sig_nhit_layer_6_cycle))[0];
    ((*sig_nhit_layer_7))[i-1] = ((*sig_nhit_layer_7_cycle))[0];
    ((*sig_nhit_layer_8))[i-1] = ((*sig_nhit_layer_8_cycle))[0];
    ((*sig_nhit_layer_9))[i-1] = ((*sig_nhit_layer_9_cycle))[0];
    ((*sig_nhit_layer_10))[i-1] = ((*sig_nhit_layer_10_cycle))[0];
    ((*sig_nhit_layer_11))[i-1] = ((*sig_nhit_layer_11_cycle))[0];
    ((*sig_nhit_layer_12))[i-1] = ((*sig_nhit_layer_12_cycle))[0];
    ((*sig_nhit_layer_13))[i-1] = ((*sig_nhit_layer_13_cycle))[0];
    ((*sig_nhit_layer_14))[i-1] = ((*sig_nhit_layer_14_cycle))[0];

    ((*mu_nhit_layer_n_0))[i-1] = ((* mu_nhit_layer_n_0_cycle))[0];
    ((*mu_nhit_layer_n_1))[i-1] = ((* mu_nhit_layer_n_1_cycle))[0];
    ((*mu_nhit_layer_n_2))[i-1] = ((* mu_nhit_layer_n_2_cycle))[0];
    ((*mu_nhit_layer_n_3))[i-1] = ((* mu_nhit_layer_n_3_cycle))[0];
    ((*mu_nhit_layer_n_4))[i-1] = ((* mu_nhit_layer_n_4_cycle))[0];
    ((*mu_nhit_layer_n_5))[i-1] = ((* mu_nhit_layer_n_5_cycle))[0];
    ((*mu_nhit_layer_n_6))[i-1] = ((* mu_nhit_layer_n_6_cycle))[0];
    ((*mu_nhit_layer_n_7))[i-1] = ((* mu_nhit_layer_n_7_cycle))[0];
    ((*mu_nhit_layer_n_8))[i-1] = ((* mu_nhit_layer_n_8_cycle))[0];
    ((*mu_nhit_layer_n_9))[i-1] = ((* mu_nhit_layer_n_9_cycle))[0];
    ((*mu_nhit_layer_n_10))[i-1] = ((* mu_nhit_layer_n_10_cycle))[0];
    ((*mu_nhit_layer_n_11))[i-1] = ((* mu_nhit_layer_n_11_cycle))[0];
    ((*mu_nhit_layer_n_12))[i-1] = ((* mu_nhit_layer_n_12_cycle))[0];
    ((*mu_nhit_layer_n_13))[i-1] = ((* mu_nhit_layer_n_13_cycle))[0];
    ((*mu_nhit_layer_n_14))[i-1] = ((* mu_nhit_layer_n_14_cycle))[0];
    ((*mu_error_nhit_layer_n_0))[i-1] = ((* mu_error_nhit_layer_n_0_cycle))[0];
    ((*mu_error_nhit_layer_n_1))[i-1] = ((* mu_error_nhit_layer_n_1_cycle))[0];
    ((*mu_error_nhit_layer_n_2))[i-1] = ((* mu_error_nhit_layer_n_2_cycle))[0];
    ((*mu_error_nhit_layer_n_3))[i-1] = ((* mu_error_nhit_layer_n_3_cycle))[0];
    ((*mu_error_nhit_layer_n_4))[i-1] = ((* mu_error_nhit_layer_n_4_cycle))[0];
    ((*mu_error_nhit_layer_n_5))[i-1] = ((* mu_error_nhit_layer_n_5_cycle))[0];
    ((*mu_error_nhit_layer_n_6))[i-1] = ((* mu_error_nhit_layer_n_6_cycle))[0];
    ((*mu_error_nhit_layer_n_7))[i-1] = ((* mu_error_nhit_layer_n_7_cycle))[0];
    ((*mu_error_nhit_layer_n_8))[i-1] = ((* mu_error_nhit_layer_n_8_cycle))[0];
    ((*mu_error_nhit_layer_n_9))[i-1] = ((* mu_error_nhit_layer_n_9_cycle))[0];
    ((*mu_error_nhit_layer_n_10))[i-1] = ((* mu_error_nhit_layer_n_10_cycle))[0];
    ((*mu_error_nhit_layer_n_11))[i-1] = ((* mu_error_nhit_layer_n_11_cycle))[0];
    ((*mu_error_nhit_layer_n_12))[i-1] = ((* mu_error_nhit_layer_n_12_cycle))[0];
    ((*mu_error_nhit_layer_n_13))[i-1] = ((* mu_error_nhit_layer_n_13_cycle))[0];
    ((*mu_error_nhit_layer_n_14))[i-1] = ((* mu_error_nhit_layer_n_14_cycle))[0];
    ((*sig_nhit_layer_n_0))[i-1] = ((* sig_nhit_layer_n_0_cycle))[0];
    ((*sig_nhit_layer_n_1))[i-1] = ((* sig_nhit_layer_n_1_cycle))[0];
    ((*sig_nhit_layer_n_2))[i-1] = ((* sig_nhit_layer_n_2_cycle))[0];
    ((*sig_nhit_layer_n_3))[i-1] = ((* sig_nhit_layer_n_3_cycle))[0];
    ((*sig_nhit_layer_n_4))[i-1] = ((* sig_nhit_layer_n_4_cycle))[0];
    ((*sig_nhit_layer_n_5))[i-1] = ((* sig_nhit_layer_n_5_cycle))[0];
    ((*sig_nhit_layer_n_6))[i-1] = ((* sig_nhit_layer_n_6_cycle))[0];
    ((*sig_nhit_layer_n_7))[i-1] = ((* sig_nhit_layer_n_7_cycle))[0];
    ((*sig_nhit_layer_n_8))[i-1] = ((* sig_nhit_layer_n_8_cycle))[0];
    ((*sig_nhit_layer_n_9))[i-1] = ((* sig_nhit_layer_n_9_cycle))[0];
    ((*sig_nhit_layer_n_10))[i-1] = ((* sig_nhit_layer_n_10_cycle))[0];
    ((*sig_nhit_layer_n_11))[i-1] = ((* sig_nhit_layer_n_11_cycle))[0];
    ((*sig_nhit_layer_n_12))[i-1] = ((* sig_nhit_layer_n_12_cycle))[0];
    ((*sig_nhit_layer_n_13))[i-1] = ((* sig_nhit_layer_n_13_cycle))[0];
    ((*sig_nhit_layer_n_14))[i-1] = ((* sig_nhit_layer_n_14_cycle))[0];
    ((*shower_nhit_max_layer))[i-1] = ((* shower_nhit_max_layer_cycle))[0];
    ((*shower_nhit_start_layer))[i-1] = ((* shower_nhit_start_layer_cycle))[0];
    ((*shower_nhit_end_layer))[i-1] = ((* shower_nhit_end_layer_cycle))[0];
    ((*shower_nhit_max))[i-1] = ((* shower_nhit_max_cycle))[0];
    ((*shower_nhit_average))[i-1] = ((* shower_nhit_average_cycle))[0];

    //weighted energy
    ((*mu_weight_layer_0))[i-1] = ((* mu_weight_layer_0_cycle))[0];
    ((*mu_weight_layer_1))[i-1] = ((* mu_weight_layer_1_cycle))[0];
    ((*mu_weight_layer_2))[i-1] = ((* mu_weight_layer_2_cycle))[0];
    ((*mu_weight_layer_3))[i-1] = ((* mu_weight_layer_3_cycle))[0];
    ((*mu_weight_layer_4))[i-1] = ((* mu_weight_layer_4_cycle))[0];
    ((*mu_weight_layer_5))[i-1] = ((* mu_weight_layer_5_cycle))[0];
    ((*mu_weight_layer_6))[i-1] = ((* mu_weight_layer_6_cycle))[0];
    ((*mu_weight_layer_7))[i-1] = ((* mu_weight_layer_7_cycle))[0];
    ((*mu_weight_layer_8))[i-1] = ((* mu_weight_layer_8_cycle))[0];
    ((*mu_weight_layer_9))[i-1] = ((* mu_weight_layer_9_cycle))[0];
    ((*mu_weight_layer_10))[i-1] = ((* mu_weight_layer_10_cycle))[0];
    ((*mu_weight_layer_11))[i-1] = ((* mu_weight_layer_11_cycle))[0];
    ((*mu_weight_layer_12))[i-1] = ((* mu_weight_layer_12_cycle))[0];
    ((*mu_weight_layer_13))[i-1] = ((* mu_weight_layer_13_cycle))[0];
    ((*mu_weight_layer_14))[i-1] = ((* mu_weight_layer_14_cycle))[0];
    ((*mu_error_weight_layer_0))[i-1] = ((* mu_error_weight_layer_0_cycle))[0];
    ((*mu_error_weight_layer_1))[i-1] = ((* mu_error_weight_layer_1_cycle))[0];
    ((*mu_error_weight_layer_2))[i-1] = ((* mu_error_weight_layer_2_cycle))[0];
    ((*mu_error_weight_layer_3))[i-1] = ((* mu_error_weight_layer_3_cycle))[0];
    ((*mu_error_weight_layer_4))[i-1] = ((* mu_error_weight_layer_4_cycle))[0];
    ((*mu_error_weight_layer_5))[i-1] = ((* mu_error_weight_layer_5_cycle))[0];
    ((*mu_error_weight_layer_6))[i-1] = ((* mu_error_weight_layer_6_cycle))[0];
    ((*mu_error_weight_layer_7))[i-1] = ((* mu_error_weight_layer_7_cycle))[0];
    ((*mu_error_weight_layer_8))[i-1] = ((* mu_error_weight_layer_8_cycle))[0];
    ((*mu_error_weight_layer_9))[i-1] = ((* mu_error_weight_layer_9_cycle))[0];
    ((*mu_error_weight_layer_10))[i-1] = ((* mu_error_weight_layer_10_cycle))[0];
    ((*mu_error_weight_layer_11))[i-1] = ((* mu_error_weight_layer_11_cycle))[0];
    ((*mu_error_weight_layer_12))[i-1] = ((* mu_error_weight_layer_12_cycle))[0];
    ((*mu_error_weight_layer_13))[i-1] = ((* mu_error_weight_layer_13_cycle))[0];
    ((*mu_error_weight_layer_14))[i-1] = ((* mu_error_weight_layer_14_cycle))[0];
    ((*sig_weight_layer_0))[i-1] = ((* sig_weight_layer_0_cycle))[0];
    ((*sig_weight_layer_1))[i-1] = ((* sig_weight_layer_1_cycle))[0];
    ((*sig_weight_layer_2))[i-1] = ((* sig_weight_layer_2_cycle))[0];
    ((*sig_weight_layer_3))[i-1] = ((* sig_weight_layer_3_cycle))[0];
    ((*sig_weight_layer_4))[i-1] = ((* sig_weight_layer_4_cycle))[0];
    ((*sig_weight_layer_5))[i-1] = ((* sig_weight_layer_5_cycle))[0];
    ((*sig_weight_layer_6))[i-1] = ((* sig_weight_layer_6_cycle))[0];
    ((*sig_weight_layer_7))[i-1] = ((* sig_weight_layer_7_cycle))[0];
    ((*sig_weight_layer_8))[i-1] = ((* sig_weight_layer_8_cycle))[0];
    ((*sig_weight_layer_9))[i-1] = ((* sig_weight_layer_9_cycle))[0];
    ((*sig_weight_layer_10))[i-1] = ((* sig_weight_layer_10_cycle))[0];
    ((*sig_weight_layer_11))[i-1] = ((* sig_weight_layer_11_cycle))[0];
    ((*sig_weight_layer_12))[i-1] = ((* sig_weight_layer_12_cycle))[0];
    ((*sig_weight_layer_13))[i-1] = ((* sig_weight_layer_13_cycle))[0];
    ((*sig_weight_layer_14))[i-1] = ((* sig_weight_layer_14_cycle))[0];

    ((*mu_weight_layer_n_0))[i-1] = ((* mu_weight_layer_n_0_cycle))[0];
    ((*mu_weight_layer_n_1))[i-1] = ((* mu_weight_layer_n_1_cycle))[0];
    ((*mu_weight_layer_n_2))[i-1] = ((* mu_weight_layer_n_2_cycle))[0];
    ((*mu_weight_layer_n_3))[i-1] = ((* mu_weight_layer_n_3_cycle))[0];
    ((*mu_weight_layer_n_4))[i-1] = ((* mu_weight_layer_n_4_cycle))[0];
    ((*mu_weight_layer_n_5))[i-1] = ((* mu_weight_layer_n_5_cycle))[0];
    ((*mu_weight_layer_n_6))[i-1] = ((* mu_weight_layer_n_6_cycle))[0];
    ((*mu_weight_layer_n_7))[i-1] = ((* mu_weight_layer_n_7_cycle))[0];
    ((*mu_weight_layer_n_8))[i-1] = ((* mu_weight_layer_n_8_cycle))[0];
    ((*mu_weight_layer_n_9))[i-1] = ((* mu_weight_layer_n_9_cycle))[0];
    ((*mu_weight_layer_n_10))[i-1] = ((* mu_weight_layer_n_10_cycle))[0];
    ((*mu_weight_layer_n_11))[i-1] = ((* mu_weight_layer_n_11_cycle))[0];
    ((*mu_weight_layer_n_12))[i-1] = ((* mu_weight_layer_n_12_cycle))[0];
    ((*mu_weight_layer_n_13))[i-1] = ((* mu_weight_layer_n_13_cycle))[0];
    ((*mu_weight_layer_n_14))[i-1] = ((* mu_weight_layer_n_14_cycle))[0];
    ((*mu_error_weight_layer_n_0))[i-1] = ((* mu_error_weight_layer_n_0_cycle))[0];
    ((*mu_error_weight_layer_n_1))[i-1] = ((* mu_error_weight_layer_n_1_cycle))[0];
    ((*mu_error_weight_layer_n_2))[i-1] = ((* mu_error_weight_layer_n_2_cycle))[0];
    ((*mu_error_weight_layer_n_3))[i-1] = ((* mu_error_weight_layer_n_3_cycle))[0];
    ((*mu_error_weight_layer_n_4))[i-1] = ((* mu_error_weight_layer_n_4_cycle))[0];
    ((*mu_error_weight_layer_n_5))[i-1] = ((* mu_error_weight_layer_n_5_cycle))[0];
    ((*mu_error_weight_layer_n_6))[i-1] = ((* mu_error_weight_layer_n_6_cycle))[0];
    ((*mu_error_weight_layer_n_7))[i-1] = ((* mu_error_weight_layer_n_7_cycle))[0];
    ((*mu_error_weight_layer_n_8))[i-1] = ((* mu_error_weight_layer_n_8_cycle))[0];
    ((*mu_error_weight_layer_n_9))[i-1] = ((* mu_error_weight_layer_n_9_cycle))[0];
    ((*mu_error_weight_layer_n_10))[i-1] = ((* mu_error_weight_layer_n_10_cycle))[0];
    ((*mu_error_weight_layer_n_11))[i-1] = ((* mu_error_weight_layer_n_11_cycle))[0];
    ((*mu_error_weight_layer_n_12))[i-1] = ((* mu_error_weight_layer_n_12_cycle))[0];
    ((*mu_error_weight_layer_n_13))[i-1] = ((* mu_error_weight_layer_n_13_cycle))[0];
    ((*mu_error_weight_layer_n_14))[i-1] = ((* mu_error_weight_layer_n_14_cycle))[0];
    ((*sig_weight_layer_n_0))[i-1] = ((* sig_weight_layer_n_0_cycle))[0];
    ((*sig_weight_layer_n_1))[i-1] = ((* sig_weight_layer_n_1_cycle))[0];
    ((*sig_weight_layer_n_2))[i-1] = ((* sig_weight_layer_n_2_cycle))[0];
    ((*sig_weight_layer_n_3))[i-1] = ((* sig_weight_layer_n_3_cycle))[0];
    ((*sig_weight_layer_n_4))[i-1] = ((* sig_weight_layer_n_4_cycle))[0];
    ((*sig_weight_layer_n_5))[i-1] = ((* sig_weight_layer_n_5_cycle))[0];
    ((*sig_weight_layer_n_6))[i-1] = ((* sig_weight_layer_n_6_cycle))[0];
    ((*sig_weight_layer_n_7))[i-1] = ((* sig_weight_layer_n_7_cycle))[0];
    ((*sig_weight_layer_n_8))[i-1] = ((* sig_weight_layer_n_8_cycle))[0];
    ((*sig_weight_layer_n_9))[i-1] = ((* sig_weight_layer_n_9_cycle))[0];
    ((*sig_weight_layer_n_10))[i-1] = ((* sig_weight_layer_n_10_cycle))[0];
    ((*sig_weight_layer_n_11))[i-1] = ((* sig_weight_layer_n_11_cycle))[0];
    ((*sig_weight_layer_n_12))[i-1] = ((* sig_weight_layer_n_12_cycle))[0];
    ((*sig_weight_layer_n_13))[i-1] = ((* sig_weight_layer_n_13_cycle))[0];
    ((*sig_weight_layer_n_14))[i-1] = ((* sig_weight_layer_n_14_cycle))[0];

    ((*shower_weight_max_layer))[i-1] = ((* shower_weight_max_layer_cycle))[0];
    ((*shower_weight_start_layer))[i-1] = ((* shower_weight_start_layer_cycle))[0];
    ((*shower_weight_end_layer))[i-1] = ((* shower_weight_end_layer_cycle))[0];
    ((*shower_weight_max))[i-1] = ((* shower_weight_max_cycle))[0];
    ((*shower_weight_average))[i-1] = ((* shower_weight_average_cycle))[0];
    }
    energies->Print();
 
    // Load all the showers
    // Here I change sig to mu_error! 
    for(int ienergy=0; ienergy < N_ENERGIES; ienergy++){
      showernhitvalues[ienergy][0]=((*mu_nhit_layer_0))[ienergy];
      showernhitvalues[ienergy][1]=((*mu_nhit_layer_1))[ienergy];
      showernhitvalues[ienergy][2]=((*mu_nhit_layer_2))[ienergy];
      showernhitvalues[ienergy][3]=((*mu_nhit_layer_3))[ienergy];
      showernhitvalues[ienergy][4]=((*mu_nhit_layer_4))[ienergy];
      showernhitvalues[ienergy][5]=((*mu_nhit_layer_5))[ienergy];
      showernhitvalues[ienergy][6]=((*mu_nhit_layer_6))[ienergy];
      showernhitvalues[ienergy][7]=((*mu_nhit_layer_7))[ienergy];
      showernhitvalues[ienergy][8]=((*mu_nhit_layer_8))[ienergy];
      showernhitvalues[ienergy][9]=((*mu_nhit_layer_9))[ienergy];
      showernhitvalues[ienergy][10]=((*mu_nhit_layer_10))[ienergy];
      showernhitvalues[ienergy][11]=((*mu_nhit_layer_11))[ienergy];
      showernhitvalues[ienergy][12]=((*mu_nhit_layer_12))[ienergy];
      showernhitvalues[ienergy][13]=((*mu_nhit_layer_13))[ienergy];
      showernhitvalues[ienergy][14]=((*mu_nhit_layer_14))[ienergy];
      showernhiterrors[ienergy][0]=((*mu_error_nhit_layer_0))[ienergy];
      showernhiterrors[ienergy][1]=((*mu_error_nhit_layer_1))[ienergy];
      showernhiterrors[ienergy][2]=((*mu_error_nhit_layer_2))[ienergy];
      showernhiterrors[ienergy][3]=((*mu_error_nhit_layer_3))[ienergy];
      showernhiterrors[ienergy][4]=((*mu_error_nhit_layer_4))[ienergy];
      showernhiterrors[ienergy][5]=((*mu_error_nhit_layer_5))[ienergy];
      showernhiterrors[ienergy][6]=((*mu_error_nhit_layer_6))[ienergy];
      showernhiterrors[ienergy][7]=((*mu_error_nhit_layer_7))[ienergy];
      showernhiterrors[ienergy][8]=((*mu_error_nhit_layer_8))[ienergy];
      showernhiterrors[ienergy][9]=((*mu_error_nhit_layer_9))[ienergy];
      showernhiterrors[ienergy][10]=((*mu_error_nhit_layer_10))[ienergy];
      showernhiterrors[ienergy][11]=((*mu_error_nhit_layer_11))[ienergy];
      showernhiterrors[ienergy][12]=((*mu_error_nhit_layer_12))[ienergy];
      showernhiterrors[ienergy][13]=((*mu_error_nhit_layer_13))[ienergy];
      showernhiterrors[ienergy][14]=((*mu_error_nhit_layer_14))[ienergy];
      cout<<"Energy ("<<((*energies))[ienergy]<<" GeV)"<<endl;
      cout<<"nhit shower"<<endl;
      cout<<((*mu_error_nhit_layer_0))[ienergy]<<" "<<((*mu_error_nhit_layer_1))[ienergy]<<" "<<((*mu_error_nhit_layer_2))[ienergy]<<" "<<((*mu_error_nhit_layer_3))[ienergy]<<" "<<
	((*mu_error_nhit_layer_4))[ienergy]<<" "<<((*mu_error_nhit_layer_5))[ienergy]<<" "<<((*mu_error_nhit_layer_6))[ienergy]<<" "<<((*mu_error_nhit_layer_7))[ienergy]<<" "<<
	((*mu_error_nhit_layer_8))[ienergy]<<" "<<((*mu_error_nhit_layer_9))[ienergy]<<" "<<((*mu_error_nhit_layer_10))[ienergy]<<" "<<((*mu_error_nhit_layer_11))[ienergy]<<" "<<
	((*mu_error_nhit_layer_12))[ienergy]<<" "<<((*mu_error_nhit_layer_13))[ienergy]<<" "<<((*mu_error_nhit_layer_14))[ienergy]<<endl;

      showernhitvalues_n[ienergy][0]=((*mu_nhit_layer_n_0))[ienergy];
      showernhitvalues_n[ienergy][1]=((*mu_nhit_layer_n_1))[ienergy];
      showernhitvalues_n[ienergy][2]=((*mu_nhit_layer_n_2))[ienergy];
      showernhitvalues_n[ienergy][3]=((*mu_nhit_layer_n_3))[ienergy];
      showernhitvalues_n[ienergy][4]=((*mu_nhit_layer_n_4))[ienergy];
      showernhitvalues_n[ienergy][5]=((*mu_nhit_layer_n_5))[ienergy];
      showernhitvalues_n[ienergy][6]=((*mu_nhit_layer_n_6))[ienergy];
      showernhitvalues_n[ienergy][7]=((*mu_nhit_layer_n_7))[ienergy];
      showernhitvalues_n[ienergy][8]=((*mu_nhit_layer_n_8))[ienergy];
      showernhitvalues_n[ienergy][9]=((*mu_nhit_layer_n_9))[ienergy];
      showernhitvalues_n[ienergy][10]=((*mu_nhit_layer_n_10))[ienergy];
      showernhitvalues_n[ienergy][11]=((*mu_nhit_layer_n_11))[ienergy];
      showernhitvalues_n[ienergy][12]=((*mu_nhit_layer_n_12))[ienergy];
      showernhitvalues_n[ienergy][13]=((*mu_nhit_layer_n_13))[ienergy];
      showernhitvalues_n[ienergy][14]=((*mu_nhit_layer_n_14))[ienergy];
      showernhiterrors_n[ienergy][0]=((*mu_error_nhit_layer_n_0))[ienergy];
      showernhiterrors_n[ienergy][1]=((*mu_error_nhit_layer_n_1))[ienergy];
      showernhiterrors_n[ienergy][2]=((*mu_error_nhit_layer_n_2))[ienergy];
      showernhiterrors_n[ienergy][3]=((*mu_error_nhit_layer_n_3))[ienergy];
      showernhiterrors_n[ienergy][4]=((*mu_error_nhit_layer_n_4))[ienergy];
      showernhiterrors_n[ienergy][5]=((*mu_error_nhit_layer_n_5))[ienergy];
      showernhiterrors_n[ienergy][6]=((*mu_error_nhit_layer_n_6))[ienergy];
      showernhiterrors_n[ienergy][7]=((*mu_error_nhit_layer_n_7))[ienergy];
      showernhiterrors_n[ienergy][8]=((*mu_error_nhit_layer_n_8))[ienergy];
      showernhiterrors_n[ienergy][9]=((*mu_error_nhit_layer_n_9))[ienergy];
      showernhiterrors_n[ienergy][10]=((*mu_error_nhit_layer_n_10))[ienergy];
      showernhiterrors_n[ienergy][11]=((*mu_error_nhit_layer_n_11))[ienergy];
      showernhiterrors_n[ienergy][12]=((*mu_error_nhit_layer_n_12))[ienergy];
      showernhiterrors_n[ienergy][13]=((*mu_error_nhit_layer_n_13))[ienergy];
      showernhiterrors_n[ienergy][14]=((*mu_error_nhit_layer_n_14))[ienergy];
      cout<<"nhit shower n"<<endl;
      cout<<((*mu_error_nhit_layer_n_0))[ienergy]<<" "<<((*mu_error_nhit_layer_n_1))[ienergy]<<" "<<((*mu_error_nhit_layer_n_2))[ienergy]<<" "<<((*mu_error_nhit_layer_n_3))[ienergy]<<" "<<
        ((*mu_error_nhit_layer_n_4))[ienergy]<<" "<<((*mu_error_nhit_layer_n_5))[ienergy]<<" "<<((*mu_error_nhit_layer_n_6))[ienergy]<<" "<<((*mu_error_nhit_layer_n_7))[ienergy]<<" "<<
        ((*mu_error_nhit_layer_n_8))[ienergy]<<" "<<((*mu_error_nhit_layer_n_9))[ienergy]<<" "<<((*mu_error_nhit_layer_n_10))[ienergy]<<" "<<((*mu_error_nhit_layer_n_11))[ienergy]<<" "<<
        ((*mu_error_nhit_layer_n_12))[ienergy]<<" "<<((*mu_error_nhit_layer_n_13))[ienergy]<<" "<<((*mu_error_nhit_layer_n_14))[ienergy]<<endl;

      showerweightvalues[ienergy][0]=((*mu_weight_layer_0))[ienergy];
      showerweightvalues[ienergy][1]=((*mu_weight_layer_1))[ienergy];
      showerweightvalues[ienergy][2]=((*mu_weight_layer_2))[ienergy];
      showerweightvalues[ienergy][3]=((*mu_weight_layer_3))[ienergy];
      showerweightvalues[ienergy][4]=((*mu_weight_layer_4))[ienergy];
      showerweightvalues[ienergy][5]=((*mu_weight_layer_5))[ienergy];
      showerweightvalues[ienergy][6]=((*mu_weight_layer_6))[ienergy];
      showerweightvalues[ienergy][7]=((*mu_weight_layer_7))[ienergy];
      showerweightvalues[ienergy][8]=((*mu_weight_layer_8))[ienergy];
      showerweightvalues[ienergy][9]=((*mu_weight_layer_9))[ienergy];
      showerweightvalues[ienergy][10]=((*mu_weight_layer_10))[ienergy];
      showerweightvalues[ienergy][11]=((*mu_weight_layer_11))[ienergy];
      showerweightvalues[ienergy][12]=((*mu_weight_layer_12))[ienergy];
      showerweightvalues[ienergy][13]=((*mu_weight_layer_13))[ienergy];
      showerweightvalues[ienergy][14]=((*mu_weight_layer_14))[ienergy];
      showerweighterrors[ienergy][0]=((*mu_error_weight_layer_0))[ienergy];
      showerweighterrors[ienergy][1]=((*mu_error_weight_layer_1))[ienergy];
      showerweighterrors[ienergy][2]=((*mu_error_weight_layer_2))[ienergy];
      showerweighterrors[ienergy][3]=((*mu_error_weight_layer_3))[ienergy];
      showerweighterrors[ienergy][4]=((*mu_error_weight_layer_4))[ienergy];
      showerweighterrors[ienergy][5]=((*mu_error_weight_layer_5))[ienergy];
      showerweighterrors[ienergy][6]=((*mu_error_weight_layer_6))[ienergy];
      showerweighterrors[ienergy][7]=((*mu_error_weight_layer_7))[ienergy];
      showerweighterrors[ienergy][8]=((*mu_error_weight_layer_8))[ienergy];
      showerweighterrors[ienergy][9]=((*mu_error_weight_layer_9))[ienergy];
      showerweighterrors[ienergy][10]=((*mu_error_weight_layer_10))[ienergy];
      showerweighterrors[ienergy][11]=((*mu_error_weight_layer_11))[ienergy];
      showerweighterrors[ienergy][12]=((*mu_error_weight_layer_12))[ienergy];
      showerweighterrors[ienergy][13]=((*mu_error_weight_layer_13))[ienergy];
      showerweighterrors[ienergy][14]=((*mu_error_weight_layer_14))[ienergy];
      cout<<"weight shower"<<endl;
      cout<<((*mu_error_weight_layer_0))[ienergy]<<" "<<((*mu_error_weight_layer_1))[ienergy]<<" "<<((*mu_error_weight_layer_2))[ienergy]<<" "<<((*mu_error_weight_layer_3))[ienergy]<<" "<<
        ((*mu_error_weight_layer_4))[ienergy]<<" "<<((*mu_error_weight_layer_5))[ienergy]<<" "<<((*mu_error_weight_layer_6))[ienergy]<<" "<<((*mu_error_weight_layer_7))[ienergy]<<" "<<
        ((*mu_error_weight_layer_8))[ienergy]<<" "<<((*mu_error_weight_layer_9))[ienergy]<<" "<<((*mu_error_weight_layer_10))[ienergy]<<" "<<((*mu_error_weight_layer_11))[ienergy]<<" "<<
        ((*mu_error_weight_layer_12))[ienergy]<<" "<<((*mu_error_weight_layer_13))[ienergy]<<" "<<((*mu_error_weight_layer_14))[ienergy]<<endl;

      showerweightvalues_n[ienergy][0]=((*mu_weight_layer_n_0))[ienergy];
      showerweightvalues_n[ienergy][1]=((*mu_weight_layer_n_1))[ienergy];
      showerweightvalues_n[ienergy][2]=((*mu_weight_layer_n_2))[ienergy];
      showerweightvalues_n[ienergy][3]=((*mu_weight_layer_n_3))[ienergy];
      showerweightvalues_n[ienergy][4]=((*mu_weight_layer_n_4))[ienergy];
      showerweightvalues_n[ienergy][5]=((*mu_weight_layer_n_5))[ienergy];
      showerweightvalues_n[ienergy][6]=((*mu_weight_layer_n_6))[ienergy];
      showerweightvalues_n[ienergy][7]=((*mu_weight_layer_n_7))[ienergy];
      showerweightvalues_n[ienergy][8]=((*mu_weight_layer_n_8))[ienergy];
      showerweightvalues_n[ienergy][9]=((*mu_weight_layer_n_9))[ienergy];
      showerweightvalues_n[ienergy][10]=((*mu_weight_layer_n_10))[ienergy];
      showerweightvalues_n[ienergy][11]=((*mu_weight_layer_n_11))[ienergy];
      showerweightvalues_n[ienergy][12]=((*mu_weight_layer_n_12))[ienergy];
      showerweightvalues_n[ienergy][13]=((*mu_weight_layer_n_13))[ienergy];
      showerweightvalues_n[ienergy][14]=((*mu_weight_layer_n_14))[ienergy];
      showerweighterrors_n[ienergy][0]=((*mu_error_weight_layer_n_0))[ienergy];
      showerweighterrors_n[ienergy][1]=((*mu_error_weight_layer_n_1))[ienergy];
      showerweighterrors_n[ienergy][2]=((*mu_error_weight_layer_n_2))[ienergy];
      showerweighterrors_n[ienergy][3]=((*mu_error_weight_layer_n_3))[ienergy];
      showerweighterrors_n[ienergy][4]=((*mu_error_weight_layer_n_4))[ienergy];
      showerweighterrors_n[ienergy][5]=((*mu_error_weight_layer_n_5))[ienergy];
      showerweighterrors_n[ienergy][6]=((*mu_error_weight_layer_n_6))[ienergy];
      showerweighterrors_n[ienergy][7]=((*mu_error_weight_layer_n_7))[ienergy];
      showerweighterrors_n[ienergy][8]=((*mu_error_weight_layer_n_8))[ienergy];
      showerweighterrors_n[ienergy][9]=((*mu_error_weight_layer_n_9))[ienergy];
      showerweighterrors_n[ienergy][10]=((*mu_error_weight_layer_n_10))[ienergy];
      showerweighterrors_n[ienergy][11]=((*mu_error_weight_layer_n_11))[ienergy];
      showerweighterrors_n[ienergy][12]=((*mu_error_weight_layer_n_12))[ienergy];
      showerweighterrors_n[ienergy][13]=((*mu_error_weight_layer_n_13))[ienergy];
      showerweighterrors_n[ienergy][14]=((*mu_error_weight_layer_n_14))[ienergy];
      cout<<"weight shower n"<<endl;
      cout<<((*mu_error_weight_layer_n_0))[ienergy]<<" "<<((*mu_error_weight_layer_n_1))[ienergy]<<" "<<((*mu_error_weight_layer_n_2))[ienergy]<<" "<<((*mu_error_weight_layer_n_3))[ienergy]<<" "<<
        ((*mu_error_weight_layer_n_4))[ienergy]<<" "<<((*mu_error_weight_layer_n_5))[ienergy]<<" "<<((*mu_error_weight_layer_n_6))[ienergy]<<" "<<((*mu_error_weight_layer_n_7))[ienergy]<<" "<<
        ((*mu_error_weight_layer_n_8))[ienergy]<<" "<<((*mu_error_weight_layer_n_9))[ienergy]<<" "<<((*mu_error_weight_layer_n_10))[ienergy]<<" "<<((*mu_error_weight_layer_n_11))[ienergy]<<" "<<
        ((*mu_error_weight_layer_n_12))[ienergy]<<" "<<((*mu_error_weight_layer_n_13))[ienergy]<<" "<<((*mu_error_weight_layer_n_14))[ienergy]<<endl;
    }

    for(int ilayer=0; ilayer < N_ECAL_LAYERS; ilayer++){
      // Index    0 1 2 3  4  5  6  7  8  9 10 11 12  13  14  15  16  17                                                         
      // Energies 2 4 6 8 10 20 30 40 50 60 70 80 90 100 125 150 175 200
      showernhitvalues_2[ilayer] = showernhitvalues[0][ilayer];
      showernhitvalues_4[ilayer] = showernhitvalues[1][ilayer];
      showernhitvalues_6[ilayer] = showernhitvalues[2][ilayer];
      showernhitvalues_10[ilayer] = showernhitvalues[4][ilayer];
      showernhitvalues_20[ilayer] = showernhitvalues[5][ilayer];
      showernhitvalues_40[ilayer] = showernhitvalues[7][ilayer];
      showernhitvalues_80[ilayer] = showernhitvalues[11][ilayer];
      showernhitvalues_150[ilayer] = showernhitvalues[15][ilayer];
      showernhiterrors_2[ilayer] = showernhiterrors[0][ilayer];
      showernhiterrors_4[ilayer] = showernhiterrors[1][ilayer];
      showernhiterrors_6[ilayer] = showernhiterrors[2][ilayer];
      showernhiterrors_10[ilayer] = showernhiterrors[4][ilayer];
      showernhiterrors_20[ilayer] = showernhiterrors[5][ilayer];
      showernhiterrors_40[ilayer] = showernhiterrors[7][ilayer];
      showernhiterrors_80[ilayer] = showernhiterrors[11][ilayer];
      showernhiterrors_150[ilayer] = showernhiterrors[15][ilayer];

      showernhitvalues_n_2[ilayer] = showernhitvalues_n[0][ilayer];
      showernhitvalues_n_4[ilayer] = showernhitvalues_n[1][ilayer];
      showernhitvalues_n_6[ilayer] = showernhitvalues_n[2][ilayer];
      showernhitvalues_n_10[ilayer] = showernhitvalues_n[4][ilayer];
      showernhitvalues_n_20[ilayer] = showernhitvalues_n[5][ilayer];
      showernhitvalues_n_40[ilayer] = showernhitvalues_n[7][ilayer];
      showernhitvalues_n_80[ilayer] = showernhitvalues_n[11][ilayer];
      showernhitvalues_n_150[ilayer] = showernhitvalues_n[15][ilayer];
      showernhiterrors_n_2[ilayer] = showernhiterrors_n[0][ilayer];
      showernhiterrors_n_4[ilayer] = showernhiterrors_n[1][ilayer];
      showernhiterrors_n_6[ilayer] = showernhiterrors_n[2][ilayer];
      showernhiterrors_n_10[ilayer] = showernhiterrors_n[4][ilayer];
      showernhiterrors_n_20[ilayer] = showernhiterrors_n[5][ilayer];
      showernhiterrors_n_40[ilayer] = showernhiterrors_n[7][ilayer];
      showernhiterrors_n_80[ilayer] = showernhiterrors_n[11][ilayer];
      showernhiterrors_n_150[ilayer] = showernhiterrors_n[15][ilayer];

      showerweightvalues_2[ilayer] = showerweightvalues[0][ilayer];
      showerweightvalues_4[ilayer] = showerweightvalues[1][ilayer];
      showerweightvalues_6[ilayer] = showerweightvalues[2][ilayer];
      showerweightvalues_10[ilayer] = showerweightvalues[4][ilayer];
      showerweightvalues_20[ilayer] = showerweightvalues[5][ilayer];
      showerweightvalues_40[ilayer] = showerweightvalues[7][ilayer];
      showerweightvalues_80[ilayer] = showerweightvalues[11][ilayer];
      showerweightvalues_150[ilayer] = showerweightvalues[15][ilayer];
      showerweighterrors_2[ilayer] = showerweighterrors[0][ilayer];
      showerweighterrors_4[ilayer] = showerweighterrors[1][ilayer];
      showerweighterrors_6[ilayer] = showerweighterrors[2][ilayer];
      showerweighterrors_10[ilayer] = showerweighterrors[4][ilayer];
      showerweighterrors_20[ilayer] = showerweighterrors[5][ilayer];
      showerweighterrors_40[ilayer] = showerweighterrors[7][ilayer];
      showerweighterrors_80[ilayer] = showerweighterrors[11][ilayer];
      showerweighterrors_150[ilayer] = showerweighterrors[15][ilayer];

      showerweightvalues_n_2[ilayer] = showerweightvalues_n[0][ilayer];
      showerweightvalues_n_4[ilayer] = showerweightvalues_n[1][ilayer];
      showerweightvalues_n_6[ilayer] = showerweightvalues_n[2][ilayer];
      showerweightvalues_n_10[ilayer] = showerweightvalues_n[4][ilayer];
      showerweightvalues_n_20[ilayer] = showerweightvalues_n[5][ilayer];
      showerweightvalues_n_40[ilayer] = showerweightvalues_n[7][ilayer];
      showerweightvalues_n_80[ilayer] = showerweightvalues_n[11][ilayer];
      showerweightvalues_n_150[ilayer] = showerweightvalues_n[15][ilayer];
      showerweighterrors_n_2[ilayer] = showerweighterrors_n[0][ilayer];
      showerweighterrors_n_4[ilayer] = showerweighterrors_n[1][ilayer];
      showerweighterrors_n_6[ilayer] = showerweighterrors_n[2][ilayer];
      showerweighterrors_n_10[ilayer] = showerweighterrors_n[4][ilayer];
      showerweighterrors_n_20[ilayer] = showerweighterrors_n[5][ilayer];
      showerweighterrors_n_40[ilayer] = showerweighterrors_n[7][ilayer];
      showerweighterrors_n_80[ilayer] = showerweighterrors_n[11][ilayer];
      showerweighterrors_n_150[ilayer] = showerweighterrors_n[15][ilayer];
    }
    
    TVectorD showernhitref_2(N_ECAL_LAYERS,showernhitvalues_2);
    shower_nhit_profile_2 = &showernhitref_2;
    TVectorD showernhitref_4(N_ECAL_LAYERS,showernhitvalues_4);
    shower_nhit_profile_4 = &showernhitref_4;
    TVectorD showernhitref_6(N_ECAL_LAYERS,showernhitvalues_6);
    shower_nhit_profile_6 = &showernhitref_6;
    TVectorD showernhitref_10(N_ECAL_LAYERS,showernhitvalues_10);
    shower_nhit_profile_10 = &showernhitref_10;
    TVectorD showernhitref_20(N_ECAL_LAYERS,showernhitvalues_20);
    shower_nhit_profile_20 = &showernhitref_20;
    TVectorD showernhitref_40(N_ECAL_LAYERS,showernhitvalues_40);
    shower_nhit_profile_40 = &showernhitref_40;
    TVectorD showernhitref_80(N_ECAL_LAYERS,showernhitvalues_80);
    shower_nhit_profile_80 = &showernhitref_80;
    TVectorD showernhitref_150(N_ECAL_LAYERS,showernhitvalues_150);
    shower_nhit_profile_150 = &showernhitref_150;
    TVectorD showernhiterrorref_2(N_ECAL_LAYERS,showernhiterrors_2);
    shower_nhit_profile_errors_2 = &showernhiterrorref_2;
    TVectorD showernhiterrorref_4(N_ECAL_LAYERS,showernhiterrors_4);
    shower_nhit_profile_errors_4 = &showernhiterrorref_4;
    TVectorD showernhiterrorref_6(N_ECAL_LAYERS,showernhiterrors_6);
    shower_nhit_profile_errors_6 = &showernhiterrorref_6;
    TVectorD showernhiterrorref_10(N_ECAL_LAYERS,showernhiterrors_10);
    shower_nhit_profile_errors_10 = &showernhiterrorref_10;
    TVectorD showernhiterrorref_20(N_ECAL_LAYERS,showernhiterrors_20);
    shower_nhit_profile_errors_20 = &showernhiterrorref_20;
    TVectorD showernhiterrorref_40(N_ECAL_LAYERS,showernhiterrors_40);
    shower_nhit_profile_errors_40 = &showernhiterrorref_40;
    TVectorD showernhiterrorref_80(N_ECAL_LAYERS,showernhiterrors_80);
    shower_nhit_profile_errors_80 = &showernhiterrorref_80;
    TVectorD showernhiterrorref_150(N_ECAL_LAYERS,showernhiterrors_150);
    shower_nhit_profile_errors_150 = &showernhiterrorref_150;

    TVectorD showernhitref_n_2(N_ECAL_LAYERS,showernhitvalues_n_2);
    shower_nhit_profile_n_2 = &showernhitref_n_2;
    TVectorD showernhitref_n_4(N_ECAL_LAYERS,showernhitvalues_n_4);
    shower_nhit_profile_n_4 = &showernhitref_n_4;
    TVectorD showernhitref_n_6(N_ECAL_LAYERS,showernhitvalues_n_6);
    shower_nhit_profile_n_6 = &showernhitref_n_6;
    TVectorD showernhitref_n_10(N_ECAL_LAYERS,showernhitvalues_n_10);
    shower_nhit_profile_n_10 = &showernhitref_n_10;
    TVectorD showernhitref_n_20(N_ECAL_LAYERS,showernhitvalues_n_20);
    shower_nhit_profile_n_20 = &showernhitref_n_20;
    TVectorD showernhitref_n_40(N_ECAL_LAYERS,showernhitvalues_n_40);
    shower_nhit_profile_n_40 = &showernhitref_n_40;
    TVectorD showernhitref_n_80(N_ECAL_LAYERS,showernhitvalues_n_80);
    shower_nhit_profile_n_80 = &showernhitref_n_80;
    TVectorD showernhitref_n_150(N_ECAL_LAYERS,showernhitvalues_n_150);
    shower_nhit_profile_n_150 = &showernhitref_n_150;
    TVectorD showernhiterrorref_n_2(N_ECAL_LAYERS,showernhiterrors_n_2);
    shower_nhit_profile_errors_n_2 = &showernhiterrorref_n_2;
    TVectorD showernhiterrorref_n_4(N_ECAL_LAYERS,showernhiterrors_n_4);
    shower_nhit_profile_errors_n_4 = &showernhiterrorref_n_4;
    TVectorD showernhiterrorref_n_6(N_ECAL_LAYERS,showernhiterrors_n_6);
    shower_nhit_profile_errors_n_6 = &showernhiterrorref_n_6;
    TVectorD showernhiterrorref_n_10(N_ECAL_LAYERS,showernhiterrors_n_10);
    shower_nhit_profile_errors_n_10 = &showernhiterrorref_n_10;
    TVectorD showernhiterrorref_n_20(N_ECAL_LAYERS,showernhiterrors_n_20);
    shower_nhit_profile_errors_n_20 = &showernhiterrorref_n_20;
    TVectorD showernhiterrorref_n_40(N_ECAL_LAYERS,showernhiterrors_n_40);
    shower_nhit_profile_errors_n_40 = &showernhiterrorref_n_40;
    TVectorD showernhiterrorref_n_80(N_ECAL_LAYERS,showernhiterrors_n_80);
    shower_nhit_profile_errors_n_80 = &showernhiterrorref_n_80;
    TVectorD showernhiterrorref_n_150(N_ECAL_LAYERS,showernhiterrors_n_150);
    shower_nhit_profile_errors_n_150 = &showernhiterrorref_n_150;

    TVectorD showerweightref_2(N_ECAL_LAYERS,showerweightvalues_2);
    shower_weight_profile_2 = &showerweightref_2;
    TVectorD showerweightref_4(N_ECAL_LAYERS,showerweightvalues_4);
    shower_weight_profile_4 = &showerweightref_4;
    TVectorD showerweightref_6(N_ECAL_LAYERS,showerweightvalues_6);
    shower_weight_profile_6 = &showerweightref_6;
    TVectorD showerweightref_10(N_ECAL_LAYERS,showerweightvalues_10);
    shower_weight_profile_10 = &showerweightref_10;
    TVectorD showerweightref_20(N_ECAL_LAYERS,showerweightvalues_20);
    shower_weight_profile_20 = &showerweightref_20;
    TVectorD showerweightref_40(N_ECAL_LAYERS,showerweightvalues_40);
    shower_weight_profile_40 = &showerweightref_40;
    TVectorD showerweightref_80(N_ECAL_LAYERS,showerweightvalues_80);
    shower_weight_profile_80 = &showerweightref_80;
    TVectorD showerweightref_150(N_ECAL_LAYERS,showerweightvalues_150);
    shower_weight_profile_150 = &showerweightref_150;
    TVectorD showerweighterrorref_2(N_ECAL_LAYERS,showerweighterrors_2);
    shower_weight_profile_errors_2 = &showerweighterrorref_2;
    TVectorD showerweighterrorref_4(N_ECAL_LAYERS,showerweighterrors_4);
    shower_weight_profile_errors_4 = &showerweighterrorref_4;
    TVectorD showerweighterrorref_6(N_ECAL_LAYERS,showerweighterrors_6);
    shower_weight_profile_errors_6 = &showerweighterrorref_6;
    TVectorD showerweighterrorref_10(N_ECAL_LAYERS,showerweighterrors_10);
    shower_weight_profile_errors_10 = &showerweighterrorref_10;
    TVectorD showerweighterrorref_20(N_ECAL_LAYERS,showerweighterrors_20);
    shower_weight_profile_errors_20 = &showerweighterrorref_20;
    TVectorD showerweighterrorref_40(N_ECAL_LAYERS,showerweighterrors_40);
    shower_weight_profile_errors_40 = &showerweighterrorref_40;
    TVectorD showerweighterrorref_80(N_ECAL_LAYERS,showerweighterrors_80);
    shower_weight_profile_errors_80 = &showerweighterrorref_80;
    TVectorD showerweighterrorref_150(N_ECAL_LAYERS,showerweighterrors_150);
    shower_weight_profile_errors_150 = &showerweighterrorref_150;

    TVectorD showerweightref_n_2(N_ECAL_LAYERS,showerweightvalues_n_2);
    shower_weight_profile_n_2 = &showerweightref_n_2;
    TVectorD showerweightref_n_4(N_ECAL_LAYERS,showerweightvalues_n_4);
    shower_weight_profile_n_4 = &showerweightref_n_4;
    TVectorD showerweightref_n_6(N_ECAL_LAYERS,showerweightvalues_n_6);
    shower_weight_profile_n_6 = &showerweightref_n_6;
    TVectorD showerweightref_n_10(N_ECAL_LAYERS,showerweightvalues_n_10);
    shower_weight_profile_n_10 = &showerweightref_n_10;
    TVectorD showerweightref_n_20(N_ECAL_LAYERS,showerweightvalues_n_20);
    shower_weight_profile_n_20 = &showerweightref_n_20;
    TVectorD showerweightref_n_40(N_ECAL_LAYERS,showerweightvalues_n_40);
    shower_weight_profile_n_40 = &showerweightref_n_40;
    TVectorD showerweightref_n_80(N_ECAL_LAYERS,showerweightvalues_n_80);
    shower_weight_profile_n_80 = &showerweightref_n_80;
    TVectorD showerweightref_n_150(N_ECAL_LAYERS,showerweightvalues_n_150);
    shower_weight_profile_n_150 = &showerweightref_n_150;
    TVectorD showerweighterrorref_n_2(N_ECAL_LAYERS,showerweighterrors_n_2);
    shower_weight_profile_errors_n_2 = &showerweighterrorref_n_2;
    TVectorD showerweighterrorref_n_4(N_ECAL_LAYERS,showerweighterrors_n_4);
    shower_weight_profile_errors_n_4 = &showerweighterrorref_n_4;
    TVectorD showerweighterrorref_n_6(N_ECAL_LAYERS,showerweighterrors_n_6);
    shower_weight_profile_errors_n_6 = &showerweighterrorref_n_6;
    TVectorD showerweighterrorref_n_10(N_ECAL_LAYERS,showerweighterrors_n_10);
    shower_weight_profile_errors_n_10 = &showerweighterrorref_n_10;
    TVectorD showerweighterrorref_n_20(N_ECAL_LAYERS,showerweighterrors_n_20);
    shower_weight_profile_errors_n_20 = &showerweighterrorref_n_20;
    TVectorD showerweighterrorref_n_40(N_ECAL_LAYERS,showerweighterrors_n_40);
    shower_weight_profile_errors_n_40 = &showerweighterrorref_n_40;
    TVectorD showerweighterrorref_n_80(N_ECAL_LAYERS,showerweighterrors_n_80);
    shower_weight_profile_errors_n_80 = &showerweighterrorref_n_80;
    TVectorD showerweighterrorref_n_150(N_ECAL_LAYERS,showerweighterrors_n_150);
    shower_weight_profile_errors_n_150 = &showerweighterrorref_n_150;

    // Prepare energy or 1/sqrt(E)
    TVectorD *energyaxis = new TVectorD(N_ENERGIES);
    double energyvalues[N_ENERGIES];
    for(int i=0; i < N_ENERGIES; i++){
      if (transformed == true) energyvalues[i]=((*energies_tr))[i];
      else energyvalues[i]=((*energies))[i];
    }
    TVectorD energyaxisref(N_ENERGIES, energyvalues);
    energyaxis = &energyaxisref;

    
    // // string e_str = to_string((int)round(energies[i_energy])) + "GeV"+icycle).c_str()));
    // // Resolution
    
    TGraph *g_res_nhit = new TGraph(*energyaxis, *res_nhit);                 
    graph_setup_add(g_res_nhit, "N hits", kBlack);
    for (int i=0;i<g_res_nhit->GetN();i++) g_res_nhit->GetY()[i] *= 100;
    TF1 *fit_res_nhit = linearfit(g_res_nhit,18,kBlack);
    
    TGraph *g_res_sume = new TGraph(*energyaxis, *res_sume);                 
    graph_setup_add(g_res_sume, "Sum energy", kRed);
    for (int i=0;i<g_res_sume->GetN();i++) g_res_sume->GetY()[i] *= 100;
    TF1 *fit_res_sume = linearfit(g_res_sume,18,kRed);
    
    TGraph *g_res_weight = new TGraph(*energyaxis, *res_weight);                 
    graph_setup_add(g_res_weight, "Weighted sum energy", kViolet);
    for (int i=0;i<g_res_weight->GetN();i++) g_res_weight->GetY()[i] *= 100;
    TF1 *fit_res_weight = linearfit(g_res_weight,18,kViolet);
    
    // Debugging tvectord access 
    for (int j = 0; j < N_ENERGIES; j++) {cout << "This energy: " << ((*energies))[j] << " GeV, 1/sqrt(E):" << ((*energies_tr))[j] << endl;}


    //Resolution plot:
    TMultiGraph *mg = new TMultiGraph();
    
    mg->Add(g_res_nhit);
    mg->Add(g_res_sume);
    mg->Add(g_res_weight);
    
    auto c = new TCanvas("c", "c", 800, 800);
    mg->Draw("AP");
    mg->SetTitle("Resolution vs energy ("+partstring+")");
    if (transformed) mg->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg->GetXaxis()->SetTitle("E[GeV]");
    mg->GetYaxis()->SetTitle("#sigma(E_{meas.})/#mu(E_{meas.}) (\%)");
    mg->GetYaxis()->SetTitleOffset(1.4);
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
    c->SaveAs("resolution_ecal_sim_"+savetrans+".png");
    }
    
    // // Linearity plots
    // We could add sim, sim_masked and data alltogether
        
    TMultiGraph *mg_lin_nhit = new TMultiGraph();
    TMultiGraph *mg_lin_sume = new TMultiGraph();
    TMultiGraph *mg_lin_weight = new TMultiGraph();
    
    TGraphErrors *g_lin_nhit = new TGraphErrors(*energyaxis, *mu_nhit, *zeros, *sig_nhit);                 
    graph_setup_add(g_lin_nhit, "n hits", kBlack);
    TF1 *fit_lin_nhit = linearfitErrors(g_lin_nhit,18,kBlack);
    
    TGraphErrors *g_lin_sume = new TGraphErrors(*energyaxis, *mu_sume, *zeros, *sig_sume);                 
    graph_setup_add(g_lin_sume, "sum energy", kRed);
    TF1 *fit_lin_sume = linearfitErrors(g_lin_sume,18,kRed);
    
    TGraphErrors *g_lin_weight = new TGraphErrors(*energyaxis, *mu_weight, *zeros, *sig_weight);                 
    graph_setup_add(g_lin_weight, "weighted sum energy", kViolet);
    TF1 *fit_lin_weight = linearfitErrors(g_lin_weight,18,kViolet);


    mg_lin_nhit->Add(g_lin_nhit);
    mg_lin_sume->Add(g_lin_sume);
    mg_lin_weight->Add(g_lin_weight);
    
    auto c_lin_nhit = new TCanvas("c_lin_nhit", "c_lin_nhit", 800, 800);
    mg_lin_nhit->Draw("AP");
    mg_lin_nhit->SetTitle("Linearity (N hits) ("+partstring+")");
    if (transformed) mg_lin_nhit->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_lin_nhit->GetXaxis()->SetTitle("E[GeV]");
    mg_lin_nhit->GetYaxis()->SetTitle("Number of hits");
    mg_lin_nhit->GetYaxis()->SetTitleOffset(1.4);
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
      c_lin_nhit->SaveAs("lin_nhit_ecal_sim_"+savetrans+".eps");
      c_lin_nhit->SaveAs("lin_nhit_ecal_sim_"+savetrans+".png");
    }

    auto c_lin_sume = new TCanvas("c_lin_sume", "c_lin_sume", 800, 800);
    mg_lin_sume->Draw("AP");
    mg_lin_sume->SetTitle("Linearity (energy) ("+partstring+")");
    if (transformed) mg_lin_sume->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_lin_sume->GetXaxis()->SetTitle("E[GeV]");
    mg_lin_sume->GetYaxis()->SetTitle("Summed energy (MIPs)");
    mg_lin_sume->GetYaxis()->SetTitleOffset(1.4);
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
      c_lin_sume->SaveAs("lin_sume_ecal_sim_"+savetrans+".png");
    }

    auto c_lin_weight = new TCanvas("c_lin_weight", "c_lin_weight", 800, 800);
    mg_lin_weight->Draw("AP");
    mg_lin_weight->SetTitle("Linearity (weighted energy) ("+partstring+")");
    if (transformed) mg_lin_weight->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_lin_weight->GetXaxis()->SetTitle("E[GeV]");
    mg_lin_weight->GetYaxis()->SetTitle("Weighted summed energy (MIPs)");
    mg_lin_weight->GetYaxis()->SetTitleOffset(1.4);
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
      c_lin_weight->SaveAs("lin_weight_ecal_sim_"+savetrans+".png");
    }

    
    
    // Moliere
    // TGraphErrors *g_mol = new TGraphErrors(N_ENERGIES, energies_tr, mol, zeros, mol_sig);   
    TGraphErrors *g_mol = new TGraphErrors(*energyaxis, *mu_mol, *zeros, *mu_error_mol);   
    graph_setup_add(g_mol, "Moliere Radius", kBlue);

    TMultiGraph *mg_mol = new TMultiGraph();
    mg_mol->Add(g_mol);

    auto c_mol = new TCanvas("c_mol", "c_mol", 800, 800);
    mg_mol->Draw("AP");
    mg_mol->SetTitle("Moliere radius vs energy ("+partstring+")");
    // mg_mol->GetXaxis()->SetTitle("Energy [GeV]");
    if (transformed) mg_mol->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_mol->GetXaxis()->SetTitle("E[GeV]");
    mg_mol->GetYaxis()->SetTitle("Radius (90%) [mm]");
    mg_mol->GetYaxis()->SetTitleOffset(1.4);
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
      c_mol->SaveAs("mol_"+savetrans+"_w_errors.eps");
      c_mol->SaveAs("mol_"+savetrans+"_w_errors.png");
    }

    // NHITS shower profile
    // Example layer 3
    TGraphErrors *g_nhit_layer_3 = new TGraphErrors(*energyaxis, *mu_nhit_layer_3, *zeros, *sig_nhit_layer_3);
    graph_setup_add(g_nhit_layer_3, "Hits in layer 3", kBlue);

    TMultiGraph *mg_nhit_layer_3 = new TMultiGraph();
    mg_nhit_layer_3->Add(g_nhit_layer_3);

    auto c_nhit_layer_3 = new TCanvas("c_nhit_layer_3", "c_nhit_layer_3", 800, 800);
    mg_nhit_layer_3->Draw("AP");
    mg_nhit_layer_3->SetTitle("Hits in layer 3 ("+partstring+")");

    if (transformed) mg_nhit_layer_3->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_nhit_layer_3->GetXaxis()->SetTitle("E[GeV]");
    mg_nhit_layer_3->GetYaxis()->SetTitle("Hits");
    mg_nhit_layer_3->GetYaxis()->SetTitleOffset(1.4);
    mg_nhit_layer_3->GetYaxis()->SetRangeUser(0,200);

    TLegend *leg_nhit_layer_3;
    if(transformed == true) leg_nhit_layer_3= new TLegend(0.45,0.65,0.65,0.85);
    else leg_nhit_layer_3= new TLegend(0.45,0.65,0.65,0.85);
    leg_nhit_layer_3->SetTextSize(0.035);
    leg_nhit_layer_3->SetTextFont(42);
    leg_nhit_layer_3->AddEntry(g_mol,"Hits in layer 3","ap");
    leg_nhit_layer_3->SetFillColor(0);
    leg_nhit_layer_3->SetLineColor(0);
    leg_nhit_layer_3->SetShadowColor(0);
    leg_nhit_layer_3->Draw();
    addCaliceLogo();
    if(save == true){
      c_nhit_layer_3->SaveAs("nhits_layer_3_"+savetrans+".eps");
      c_nhit_layer_3->SaveAs("nhits_layer_3_"+savetrans+".png");
    }

    // Shower different energies
    TGraphErrors *g_shower_nhit_2 = new TGraphErrors(*layeraxis, *shower_nhit_profile_2 , *zeros, *shower_nhit_profile_errors_2);
    graph_setup_add(g_shower_nhit_2, "Shower profile (2 GeV)", kBlack);
    TGraphErrors *g_shower_nhit_4 = new TGraphErrors(*layeraxis, *shower_nhit_profile_4 , *zeros, *shower_nhit_profile_errors_4);
    graph_setup_add(g_shower_nhit_4, "Shower profile (4 GeV)", kCyan+3);
    TGraphErrors *g_shower_nhit_6 = new TGraphErrors(*layeraxis, *shower_nhit_profile_6 , *zeros, *shower_nhit_profile_errors_6);
    graph_setup_add(g_shower_nhit_6, "Shower profile (6 GeV)", kCyan-6);
    TGraphErrors *g_shower_nhit_10 = new TGraphErrors(*layeraxis, *shower_nhit_profile_10 , *zeros, *shower_nhit_profile_errors_10);
    graph_setup_add(g_shower_nhit_10, "Shower profile (10 GeV)", kCyan-7);
    TGraphErrors *g_shower_nhit_20 = new TGraphErrors(*layeraxis, *shower_nhit_profile_20 , *zeros, *shower_nhit_profile_errors_20);
    graph_setup_add(g_shower_nhit_20, "Shower profile (20 GeV)", kRed-9);
    TGraphErrors *g_shower_nhit_40 = new TGraphErrors(*layeraxis, *shower_nhit_profile_40 , *zeros, *shower_nhit_profile_errors_40);
    graph_setup_add(g_shower_nhit_40, "Shower profile (40 GeV)", kRed-7);
    TGraphErrors *g_shower_nhit_80 = new TGraphErrors(*layeraxis, *shower_nhit_profile_80 , *zeros, *shower_nhit_profile_errors_80);
    graph_setup_add(g_shower_nhit_80, "Shower profile (80 GeV)", kRed-3);
    TGraphErrors *g_shower_nhit_150 = new TGraphErrors(*layeraxis, *shower_nhit_profile_150 , *zeros, *shower_nhit_profile_errors_150);
    graph_setup_add(g_shower_nhit_150, "Shower profile (150 GeV)", kRed+2);

    TMultiGraph *mg_shower_nhit_energies = new TMultiGraph();
    mg_shower_nhit_energies->Add(g_shower_nhit_2);
    mg_shower_nhit_energies->Add(g_shower_nhit_4);
    mg_shower_nhit_energies->Add(g_shower_nhit_6);
    mg_shower_nhit_energies->Add(g_shower_nhit_10);
    mg_shower_nhit_energies->Add(g_shower_nhit_20);
    mg_shower_nhit_energies->Add(g_shower_nhit_40);
    mg_shower_nhit_energies->Add(g_shower_nhit_80);
    mg_shower_nhit_energies->Add(g_shower_nhit_150);

    auto c_shower_nhit_energies = new TCanvas("c_shower_nhit_energies", "c_shower_nhit_energies", 800, 800);
    mg_shower_nhit_energies->Draw("ALP");
    mg_shower_nhit_energies->SetTitle("Shower profile ("+partstring+")");

    if (transformed) mg_shower_nhit_energies->GetXaxis()->SetTitle("Layer");
    else mg_shower_nhit_energies->GetXaxis()->SetTitle("Layer");
    mg_shower_nhit_energies->GetYaxis()->SetTitle("Hits");
    mg_shower_nhit_energies->GetYaxis()->SetTitleOffset(1.4);
    mg_shower_nhit_energies->GetYaxis()->SetRangeUser(0,150);

    TLegend *leg_shower_nhit_energies;
    if(transformed == true) leg_shower_nhit_energies= new TLegend(0.12,0.60,0.32,0.85);
    else leg_shower_nhit_energies= new TLegend(0.12,0.60,0.32,0.85);
    leg_shower_nhit_energies->SetTextSize(0.035);
    leg_shower_nhit_energies->SetTextFont(42);
    leg_shower_nhit_energies->AddEntry(g_shower_nhit_2,"2 GeV","ap");
    leg_shower_nhit_energies->AddEntry(g_shower_nhit_4,"4 GeV","ap");
    leg_shower_nhit_energies->AddEntry(g_shower_nhit_6,"6 GeV","ap");
    leg_shower_nhit_energies->AddEntry(g_shower_nhit_10,"10 GeV","ap");
    leg_shower_nhit_energies->AddEntry(g_shower_nhit_20,"20 GeV","ap");
    leg_shower_nhit_energies->AddEntry(g_shower_nhit_40,"40 GeV","ap");
    leg_shower_nhit_energies->AddEntry(g_shower_nhit_80,"80 GeV","ap");
    leg_shower_nhit_energies->AddEntry(g_shower_nhit_150,"150 GeV","ap");
    leg_shower_nhit_energies->SetFillColor(0);
    leg_shower_nhit_energies->SetLineColor(0);
    leg_shower_nhit_energies->SetShadowColor(0);
    leg_shower_nhit_energies->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_nhit_energies->SaveAs("shower_nhit_energies_"+savetrans+"_w_errors.eps");
      c_shower_nhit_energies->SaveAs("shower_nhit_energies_"+savetrans+"_w_errors.png");
    }


    // Shower profile normalized
    // Example layer 3                                                                                                                                          
    TGraphErrors *g_nhit_layer_n_3 = new TGraphErrors(*energyaxis, *mu_nhit_layer_n_3, *zeros, *sig_nhit_layer_n_3);
    graph_setup_add(g_nhit_layer_n_3, "Hits in layer 3 (normalized)", kBlue);

    TMultiGraph *mg_nhit_layer_n_3 = new TMultiGraph();
    mg_nhit_layer_n_3->Add(g_nhit_layer_n_3);

    auto c_nhit_layer_n_3 = new TCanvas("c_nhit_layer_n_3", "c_nhit_layer_n_3", 800, 800);
    mg_nhit_layer_n_3->Draw("AP");
    mg_nhit_layer_n_3->SetTitle("Hits in layer 3 (normalized) ("+partstring+")");

    if (transformed) mg_nhit_layer_n_3->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_nhit_layer_n_3->GetXaxis()->SetTitle("E[GeV]");
    mg_nhit_layer_n_3->GetYaxis()->SetTitle("Hits");
    mg_nhit_layer_n_3->GetYaxis()->SetTitleOffset(1.4);
    mg_nhit_layer_n_3->GetYaxis()->SetRangeUser(0.,0.3);

    TLegend *leg_nhit_layer_n_3;
    if(transformed == true) leg_nhit_layer_n_3= new TLegend(0.45,0.65,0.65,0.85);
    else leg_nhit_layer_n_3= new TLegend(0.45,0.65,0.65,0.85);
    leg_nhit_layer_n_3->SetTextSize(0.035);
    leg_nhit_layer_n_3->SetTextFont(42);
    leg_nhit_layer_n_3->AddEntry(g_mol,"Hits in layer 3 (normalized)","ap");
    leg_nhit_layer_n_3->SetFillColor(0);
    leg_nhit_layer_n_3->SetLineColor(0);
    leg_nhit_layer_n_3->SetShadowColor(0);
    leg_nhit_layer_n_3->Draw();
    addCaliceLogo();
    if(save == true){
      c_nhit_layer_n_3->SaveAs("nhit_layer_n_3_"+savetrans+"_w_errors.eps");
      c_nhit_layer_n_3->SaveAs("nhit_layer_n_3_"+savetrans+"_w_errors.png");
    }

    // Shower nhit normalized for different energies                                                                                                                                                           
    TGraphErrors *g_shower_nhit_n_2 = new TGraphErrors(*layeraxis, *shower_nhit_profile_n_2 , *zeros, *shower_nhit_profile_errors_n_2);
    graph_setup_add(g_shower_nhit_n_2, "Shower profile w (2 GeV)", kBlack);
    TGraphErrors *g_shower_nhit_n_4 = new TGraphErrors(*layeraxis, *shower_nhit_profile_n_4 , *zeros, *shower_nhit_profile_errors_n_4);
    graph_setup_add(g_shower_nhit_n_4, "Shower profile w (4 GeV)", kCyan+3);
    TGraphErrors *g_shower_nhit_n_6 = new TGraphErrors(*layeraxis, *shower_nhit_profile_n_6 , *zeros, *shower_nhit_profile_errors_n_6);
    graph_setup_add(g_shower_nhit_n_6, "Shower profile w (6 GeV)", kCyan-6);
    TGraphErrors *g_shower_nhit_n_10 = new TGraphErrors(*layeraxis, *shower_nhit_profile_n_10 , *zeros, *shower_nhit_profile_errors_n_10);
    graph_setup_add(g_shower_nhit_n_10, "Shower profile w (10 GeV)", kCyan-7);
    TGraphErrors *g_shower_nhit_n_20 = new TGraphErrors(*layeraxis, *shower_nhit_profile_n_20 , *zeros, *shower_nhit_profile_errors_n_20);
    graph_setup_add(g_shower_nhit_n_20, "Shower profile w (20 GeV)", kRed-9);
    TGraphErrors *g_shower_nhit_n_40 = new TGraphErrors(*layeraxis, *shower_nhit_profile_n_40 , *zeros, *shower_nhit_profile_errors_n_40);
    graph_setup_add(g_shower_nhit_n_40, "Shower profile w (40 GeV)", kRed-7);
    TGraphErrors *g_shower_nhit_n_80 = new TGraphErrors(*layeraxis, *shower_nhit_profile_n_80 , *zeros, *shower_nhit_profile_errors_n_80);
    graph_setup_add(g_shower_nhit_n_80, "Shower profile w (80 GeV)", kRed-3);
    TGraphErrors *g_shower_nhit_n_150 = new TGraphErrors(*layeraxis, *shower_nhit_profile_n_150 , *zeros, *shower_nhit_profile_errors_n_150);
    graph_setup_add(g_shower_nhit_n_150, "Shower profile w (150 GeV)", kRed+2);

    TMultiGraph *mg_shower_nhit_energies_n = new TMultiGraph();
    mg_shower_nhit_energies_n->Add(g_shower_nhit_n_2);
    mg_shower_nhit_energies_n->Add(g_shower_nhit_n_4);
    mg_shower_nhit_energies_n->Add(g_shower_nhit_n_6);
    mg_shower_nhit_energies_n->Add(g_shower_nhit_n_10);
    mg_shower_nhit_energies_n->Add(g_shower_nhit_n_20);
    mg_shower_nhit_energies_n->Add(g_shower_nhit_n_40);
    mg_shower_nhit_energies_n->Add(g_shower_nhit_n_80);
    mg_shower_nhit_energies_n->Add(g_shower_nhit_n_150);

    auto c_shower_nhit_energies_n = new TCanvas("c_shower_nhit_energies_n", "c_shower_nhit_energies_n", 800, 800);
    mg_shower_nhit_energies_n->Draw("ALP");
    mg_shower_nhit_energies_n->SetTitle("Shower profile (normalized) ("+partstring+")");

    if (transformed) mg_shower_nhit_energies_n->GetXaxis()->SetTitle("Layer");
    else mg_shower_nhit_energies_n->GetXaxis()->SetTitle("Layer");
    mg_shower_nhit_energies_n->GetYaxis()->SetTitle("Hits (normalized)");
    mg_shower_nhit_energies_n->GetYaxis()->SetTitleOffset(1.4);
    mg_shower_nhit_energies_n->GetYaxis()->SetRangeUser(0.,0.3);

    TLegend *leg_shower_nhit_energies_n;
    if(transformed == true) leg_shower_nhit_energies_n= new TLegend(0.12,0.60,0.32,0.85);
    else leg_shower_nhit_energies_n= new TLegend(0.12,0.60,0.32,0.85);
    leg_shower_nhit_energies_n->SetTextSize(0.035);
    leg_shower_nhit_energies_n->SetTextFont(42);
    leg_shower_nhit_energies_n->AddEntry(g_shower_nhit_n_2,"2 GeV","ap");
    leg_shower_nhit_energies_n->AddEntry(g_shower_nhit_n_4,"4 GeV","ap");
    leg_shower_nhit_energies_n->AddEntry(g_shower_nhit_n_6,"6 GeV","ap");
    leg_shower_nhit_energies_n->AddEntry(g_shower_nhit_n_10,"10 GeV","ap");
    leg_shower_nhit_energies_n->AddEntry(g_shower_nhit_n_20,"20 GeV","ap");
    leg_shower_nhit_energies_n->AddEntry(g_shower_nhit_n_40,"40 GeV","ap");
    leg_shower_nhit_energies_n->AddEntry(g_shower_nhit_n_80,"80 GeV","ap");
    leg_shower_nhit_energies_n->AddEntry(g_shower_nhit_n_150,"150 GeV","ap");
    leg_shower_nhit_energies_n->SetFillColor(0);
    leg_shower_nhit_energies_n->SetLineColor(0);
    leg_shower_nhit_energies_n->SetShadowColor(0);
    leg_shower_nhit_energies_n->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_nhit_energies_n->SaveAs("shower_nhit_energies_n_"+savetrans+"_w_errors.eps");
      c_shower_nhit_energies_n->SaveAs("shower_nhit_energies_n_"+savetrans+"_w_errors.png");
    }

    // Shower max, start and end
    // Shower nhit max
    TGraphErrors *g_shower_nhit_max = new TGraphErrors(*energyaxis, *shower_nhit_max, *zeros, *zeros);
    graph_setup_add(g_shower_nhit_max, "Shower maximum (n hits)", kBlue);

    TMultiGraph *mg_shower_nhit_max = new TMultiGraph();
    mg_shower_nhit_max->Add(g_shower_nhit_max);

    auto c_shower_nhit_max = new TCanvas("c_shower_nhit_max", "c_shower_nhit_max", 800, 800);
    mg_shower_nhit_max->Draw("AP");
    mg_shower_nhit_max->SetTitle("Shower max (n hits) ("+partstring+")");

    if (transformed) mg_shower_nhit_max->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_nhit_max->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_nhit_max->GetYaxis()->SetTitle("Layer");
    mg_shower_nhit_max->GetYaxis()->SetTitleOffset(1.4);
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
      c_shower_nhit_max->SaveAs("shower_nhit_max_"+savetrans+".eps");
      c_shower_nhit_max->SaveAs("shower_nhit_max_"+savetrans+".png");
    }

    // Shower nhit average
    TGraphErrors *g_shower_nhit_average = new TGraphErrors(*energyaxis, *shower_nhit_average, *zeros, *zeros);
    graph_setup_add(g_shower_nhit_average, "Shower average (n hits)", kBlue);

    TMultiGraph *mg_shower_nhit_average = new TMultiGraph();
    mg_shower_nhit_average->Add(g_shower_nhit_average);

    auto c_shower_nhit_average = new TCanvas("c_shower_nhit_average", "c_shower_nhit_average", 800, 800);
    mg_shower_nhit_average->Draw("AP");
    mg_shower_nhit_average->SetTitle("Shower average (n hits) ("+partstring+")");

    if (transformed) mg_shower_nhit_average->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_nhit_average->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_nhit_average->GetYaxis()->SetTitle("Layer");
    mg_shower_nhit_average->GetYaxis()->SetTitleOffset(1.4);
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
      c_shower_nhit_average->SaveAs("shower_nhit_average_"+savetrans+".eps");
      c_shower_nhit_average->SaveAs("shower_nhit_average_"+savetrans+".png");
    }

    // Shower nhit max layer
    TGraphErrors *g_shower_nhit_max_layer = new TGraphErrors(*energyaxis, *shower_nhit_max_layer, *zeros, *zeros);
    graph_setup_add(g_shower_nhit_max_layer, "Shower maximum (n hits)", kBlue);

    TMultiGraph *mg_shower_nhit_max_layer = new TMultiGraph();
    mg_shower_nhit_max_layer->Add(g_shower_nhit_max_layer);

    auto c_shower_nhit_max_layer = new TCanvas("c_shower_nhit_max_layer", "c_shower_nhit_max_layer", 800, 800);
    mg_shower_nhit_max_layer->Draw("AP");
    mg_shower_nhit_max_layer->SetTitle("N hits shower max layer (avg.) ("+partstring+")");

    if (transformed) mg_shower_nhit_max_layer->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_nhit_max_layer->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_nhit_max_layer->GetYaxis()->SetTitle("Layer");
    mg_shower_nhit_max_layer->GetYaxis()->SetTitleOffset(1.4);
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
      c_shower_nhit_max_layer->SaveAs("shower_nhit_max_layer_"+savetrans+".eps");
      c_shower_nhit_max_layer->SaveAs("shower_nhit_max_layer_"+savetrans+".png");
    }

    // Shower nhit start layer
    TGraphErrors *g_shower_nhit_start_layer = new TGraphErrors(*energyaxis, *shower_nhit_start_layer, *zeros, *zeros);
    graph_setup_add(g_shower_nhit_start_layer, "Shower startimum (n hits)", kBlue);

    TMultiGraph *mg_shower_nhit_start_layer = new TMultiGraph();
    mg_shower_nhit_start_layer->Add(g_shower_nhit_start_layer);

    auto c_shower_nhit_start_layer = new TCanvas("c_shower_nhit_start_layer", "c_shower_nhit_start_layer", 800, 800);
    mg_shower_nhit_start_layer->Draw("AP");
    mg_shower_nhit_start_layer->SetTitle("N hits shower start layer (avg.) ("+partstring+")");

    if (transformed) mg_shower_nhit_start_layer->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_nhit_start_layer->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_nhit_start_layer->GetYaxis()->SetTitle("Layer");
    mg_shower_nhit_start_layer->GetYaxis()->SetTitleOffset(1.4);
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
      c_shower_nhit_start_layer->SaveAs("shower_nhit_start_layer_"+savetrans+".eps");
      c_shower_nhit_start_layer->SaveAs("shower_nhit_start_layer_"+savetrans+".png");
    }

    // Shower nhit end layer                                                                                                                                                                                                            
    TGraphErrors *g_shower_nhit_end_layer = new TGraphErrors(*energyaxis, *shower_nhit_end_layer, *zeros, *zeros);
    graph_setup_add(g_shower_nhit_end_layer, "Shower endimum (n hits)", kBlue);

    TMultiGraph *mg_shower_nhit_end_layer = new TMultiGraph();
    mg_shower_nhit_end_layer->Add(g_shower_nhit_end_layer);

    auto c_shower_nhit_end_layer = new TCanvas("c_shower_nhit_end_layer", "c_shower_nhit_end_layer", 800, 800);
    mg_shower_nhit_end_layer->Draw("AP");
    mg_shower_nhit_end_layer->SetTitle("N hits shower end layer (avg.) ("+partstring+")");

    if (transformed) mg_shower_nhit_end_layer->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_nhit_end_layer->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_nhit_end_layer->GetYaxis()->SetTitle("Layer");
    mg_shower_nhit_end_layer->GetYaxis()->SetTitleOffset(1.4);
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
      c_shower_nhit_end_layer->SaveAs("shower_nhit_end_layer_"+savetrans+".eps");
      c_shower_nhit_end_layer->SaveAs("shower_nhit_end_layer_"+savetrans+".png");
    }

    //___________---------_______________________-------______--_-_-____-__--__-_-_----_-_-_--_-__---


    // WEIGHT shower profile
    // Example layer 3
    TGraphErrors *g_weight_layer_3 = new TGraphErrors(*energyaxis, *mu_weight_layer_3, *zeros, *sig_weight_layer_3);
    graph_setup_add(g_weight_layer_3, "Weighted energy in layer 3", kBlue);

    TMultiGraph *mg_weight_layer_3 = new TMultiGraph();
    mg_weight_layer_3->Add(g_weight_layer_3);

    auto c_weight_layer_3 = new TCanvas("c_weight_layer_3", "c_weight_layer_3", 800, 800);
    mg_weight_layer_3->Draw("AP");
    mg_weight_layer_3->SetTitle("Weighted energy in layer 3 ("+partstring+")");

    if (transformed) mg_weight_layer_3->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_weight_layer_3->GetXaxis()->SetTitle("E[GeV]");
    mg_weight_layer_3->GetYaxis()->SetTitle("Weighted energy (MIPs)");
    mg_weight_layer_3->GetYaxis()->SetTitleOffset(1.4);
    mg_weight_layer_3->GetYaxis()->SetRangeUser(0,200);

    TLegend *leg_weight_layer_3;
    if(transformed == true) leg_weight_layer_3= new TLegend(0.45,0.65,0.65,0.85);
    else leg_weight_layer_3= new TLegend(0.45,0.65,0.65,0.85);
    leg_weight_layer_3->SetTextSize(0.035);
    leg_weight_layer_3->SetTextFont(42);
    leg_weight_layer_3->AddEntry(g_mol,"Weighted energy in layer 3","ap");
    leg_weight_layer_3->SetFillColor(0);
    leg_weight_layer_3->SetLineColor(0);
    leg_weight_layer_3->SetShadowColor(0);
    leg_weight_layer_3->Draw();
    addCaliceLogo();
    if(save == true){
      c_weight_layer_3->SaveAs("weight_layer_3_"+savetrans+"_w_errors.eps");
      c_weight_layer_3->SaveAs("weight_layer_3_"+savetrans+"_w_errors.png");
    }

    // Shower different energies
    TGraphErrors *g_shower_weight_2 = new TGraphErrors(*layeraxis, *shower_weight_profile_2 , *zeros, *shower_weight_profile_errors_2);
    graph_setup_add(g_shower_weight_2, "Shower profile (2 GeV)", kBlack);
    TGraphErrors *g_shower_weight_4 = new TGraphErrors(*layeraxis, *shower_weight_profile_4 , *zeros, *shower_weight_profile_errors_4);
    graph_setup_add(g_shower_weight_4, "Shower profile (4 GeV)", kCyan+3);
    TGraphErrors *g_shower_weight_6 = new TGraphErrors(*layeraxis, *shower_weight_profile_6 , *zeros, *shower_weight_profile_errors_6);
    graph_setup_add(g_shower_weight_6, "Shower profile (6 GeV)", kCyan-6);
    TGraphErrors *g_shower_weight_10 = new TGraphErrors(*layeraxis, *shower_weight_profile_10 , *zeros, *shower_weight_profile_errors_10);
    graph_setup_add(g_shower_weight_10, "Shower profile (10 GeV)", kCyan-7);
    TGraphErrors *g_shower_weight_20 = new TGraphErrors(*layeraxis, *shower_weight_profile_20 , *zeros, *shower_weight_profile_errors_20);
    graph_setup_add(g_shower_weight_20, "Shower profile (20 GeV)", kRed-9);
    TGraphErrors *g_shower_weight_40 = new TGraphErrors(*layeraxis, *shower_weight_profile_40 , *zeros, *shower_weight_profile_errors_40);
    graph_setup_add(g_shower_weight_40, "Shower profile (40 GeV)", kRed-7);
    TGraphErrors *g_shower_weight_80 = new TGraphErrors(*layeraxis, *shower_weight_profile_80 , *zeros, *shower_weight_profile_errors_80);
    graph_setup_add(g_shower_weight_80, "Shower profile (80 GeV)", kRed-3);
    TGraphErrors *g_shower_weight_150 = new TGraphErrors(*layeraxis, *shower_weight_profile_150 , *zeros, *shower_weight_profile_errors_150);
    graph_setup_add(g_shower_weight_150, "Shower profile (150 GeV)", kRed+2);

    TMultiGraph *mg_shower_weight_energies = new TMultiGraph();
    mg_shower_weight_energies->Add(g_shower_weight_2);
    mg_shower_weight_energies->Add(g_shower_weight_4);
    mg_shower_weight_energies->Add(g_shower_weight_6);
    mg_shower_weight_energies->Add(g_shower_weight_10);
    mg_shower_weight_energies->Add(g_shower_weight_20);
    mg_shower_weight_energies->Add(g_shower_weight_40);
    mg_shower_weight_energies->Add(g_shower_weight_80);
    mg_shower_weight_energies->Add(g_shower_weight_150);

    auto c_shower_weight_energies = new TCanvas("c_shower_weight_energies", "c_shower_weight_energies", 800, 800);
    mg_shower_weight_energies->Draw("ALP");
    mg_shower_weight_energies->SetTitle("Shower profile ("+partstring+")");
    mg_shower_weight_energies->GetYaxis()->SetMaxDigits(3);

    if (transformed) mg_shower_weight_energies->GetXaxis()->SetTitle("Layer");
    else mg_shower_weight_energies->GetXaxis()->SetTitle("Layer");
    mg_shower_weight_energies->GetYaxis()->SetTitle("Weighted energy (MIPs)");
    mg_shower_weight_energies->GetYaxis()->SetTitleOffset(1.4);
    mg_shower_weight_energies->GetYaxis()->SetRangeUser(0,7000);

    TLegend *leg_shower_weight_energies;
    if(transformed == true) leg_shower_weight_energies= new TLegend(0.12,0.60,0.32,0.85);
    else leg_shower_weight_energies= new TLegend(0.12,0.60,0.32,0.85);
    leg_shower_weight_energies->SetTextSize(0.035);
    leg_shower_weight_energies->SetTextFont(42);
    leg_shower_weight_energies->AddEntry(g_shower_weight_2,"2 GeV","ap");
    leg_shower_weight_energies->AddEntry(g_shower_weight_4,"4 GeV","ap");
    leg_shower_weight_energies->AddEntry(g_shower_weight_6,"6 GeV","ap");
    leg_shower_weight_energies->AddEntry(g_shower_weight_10,"10 GeV","ap");
    leg_shower_weight_energies->AddEntry(g_shower_weight_20,"20 GeV","ap");
    leg_shower_weight_energies->AddEntry(g_shower_weight_40,"40 GeV","ap");
    leg_shower_weight_energies->AddEntry(g_shower_weight_80,"80 GeV","ap");
    leg_shower_weight_energies->AddEntry(g_shower_weight_150,"150 GeV","ap");
    leg_shower_weight_energies->SetFillColor(0);
    leg_shower_weight_energies->SetLineColor(0);
    leg_shower_weight_energies->SetShadowColor(0);
    leg_shower_weight_energies->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_weight_energies->SaveAs("shower_weight_energies_"+savetrans+"_w_errors.eps");
      c_shower_weight_energies->SaveAs("shower_weight_energies_"+savetrans+"_w_errors.png");
    }


    // Shower profile normalized
    // Example layer 3
    TGraphErrors *g_weights_layer_n_3 = new TGraphErrors(*energyaxis, *mu_weight_layer_n_3, *zeros, *zeros);
    graph_setup_add(g_weights_layer_n_3, "Weighted energy in layer 3 (normalized)", kBlue);

    TMultiGraph *mg_weights_layer_n_3 = new TMultiGraph();
    mg_weights_layer_n_3->Add(g_weights_layer_n_3);

    auto c_weight_layer_n_3 = new TCanvas("c_weight_layer_n_3", "c_weight_layer_n_3", 800, 800);
    mg_weights_layer_n_3->Draw("AP");
    mg_weights_layer_n_3->SetTitle("W. energy in layer 3 (normalized) ("+partstring+")");

    if (transformed) mg_weights_layer_n_3->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_weights_layer_n_3->GetXaxis()->SetTitle("E[GeV]");
    mg_weights_layer_n_3->GetYaxis()->SetTitle("Weighted energy (MIPs)");
    mg_weights_layer_n_3->GetYaxis()->SetTitleOffset(1.4);
    mg_weights_layer_n_3->GetYaxis()->SetRangeUser(0.,0.3);

    TLegend *leg_weight_layer_n_3;
    if(transformed == true) leg_weight_layer_n_3= new TLegend(0.45,0.65,0.65,0.85);
    else leg_weight_layer_n_3= new TLegend(0.45,0.65,0.65,0.85);
    leg_weight_layer_n_3->SetTextSize(0.035);
    leg_weight_layer_n_3->SetTextFont(42);
    leg_weight_layer_n_3->AddEntry(g_mol,"Weighted energy in layer 3 (normalized)","ap");
    leg_weight_layer_n_3->SetFillColor(0);
    leg_weight_layer_n_3->SetLineColor(0);
    leg_weight_layer_n_3->SetShadowColor(0);
    leg_weight_layer_n_3->Draw();
    addCaliceLogo();
    if(save == true){
      c_weight_layer_n_3->SaveAs("weight_layer_n_3_"+savetrans+"_w_errors.eps");
      c_weight_layer_n_3->SaveAs("weight_layer_n_3_"+savetrans+"_w_errors.png");
    }

    // Shower weight normalized for different energies
    TGraphErrors *g_shower_weight_n_2 = new TGraphErrors(*layeraxis, *shower_weight_profile_n_2 , *zeros, *shower_weight_profile_errors_n_2);
    graph_setup_add(g_shower_weight_n_2, "Shower profile w (2 GeV)", kBlack);
    TGraphErrors *g_shower_weight_n_4 = new TGraphErrors(*layeraxis, *shower_weight_profile_n_4 , *zeros, *shower_weight_profile_errors_n_4);
    graph_setup_add(g_shower_weight_n_4, "Shower profile w (4 GeV)", kCyan+3);
    TGraphErrors *g_shower_weight_n_6 = new TGraphErrors(*layeraxis, *shower_weight_profile_n_6 , *zeros, *shower_weight_profile_errors_n_6);
    graph_setup_add(g_shower_weight_n_6, "Shower profile w (6 GeV)", kCyan-6);
    TGraphErrors *g_shower_weight_n_10 = new TGraphErrors(*layeraxis, *shower_weight_profile_n_10 , *zeros, *shower_weight_profile_errors_n_10);
    graph_setup_add(g_shower_weight_n_10, "Shower profile w (10 GeV)", kCyan-7);
    TGraphErrors *g_shower_weight_n_20 = new TGraphErrors(*layeraxis, *shower_weight_profile_n_20 , *zeros, *shower_weight_profile_errors_n_20);
    graph_setup_add(g_shower_weight_n_20, "Shower profile w (20 GeV)", kRed-9);
    TGraphErrors *g_shower_weight_n_40 = new TGraphErrors(*layeraxis, *shower_weight_profile_n_40 , *zeros, *shower_weight_profile_errors_n_40);
    graph_setup_add(g_shower_weight_n_40, "Shower profile w (40 GeV)", kRed-7);
    TGraphErrors *g_shower_weight_n_80 = new TGraphErrors(*layeraxis, *shower_weight_profile_n_80 , *zeros, *shower_weight_profile_errors_n_80);
    graph_setup_add(g_shower_weight_n_80, "Shower profile w (80 GeV)", kRed-3);
    TGraphErrors *g_shower_weight_n_150 = new TGraphErrors(*layeraxis, *shower_weight_profile_n_150 , *zeros, *shower_weight_profile_errors_n_150);
    graph_setup_add(g_shower_weight_n_150, "Shower profile w (150 GeV)", kRed+2);

    TMultiGraph *mg_shower_weight_energies_n = new TMultiGraph();
    mg_shower_weight_energies_n->Add(g_shower_weight_n_2);
    mg_shower_weight_energies_n->Add(g_shower_weight_n_4);
    mg_shower_weight_energies_n->Add(g_shower_weight_n_6);
    mg_shower_weight_energies_n->Add(g_shower_weight_n_10);
    mg_shower_weight_energies_n->Add(g_shower_weight_n_20);
    mg_shower_weight_energies_n->Add(g_shower_weight_n_40);
    mg_shower_weight_energies_n->Add(g_shower_weight_n_80);
    mg_shower_weight_energies_n->Add(g_shower_weight_n_150);

    auto c_shower_weight_energies_n = new TCanvas("c_shower_weight_energies_n", "c_shower_weight_energies_n", 800, 800);
    mg_shower_weight_energies_n->Draw("ALP");
    mg_shower_weight_energies_n->SetTitle("Shower profile (normalized) ("+partstring+")");

    if (transformed) mg_shower_weight_energies_n->GetXaxis()->SetTitle("Layer");
    else mg_shower_weight_energies_n->GetXaxis()->SetTitle("Layer");
    mg_shower_weight_energies_n->GetYaxis()->SetTitle("Weighted energy (normalized)");
    mg_shower_weight_energies_n->GetYaxis()->SetTitleOffset(1.4);
    mg_shower_weight_energies_n->GetYaxis()->SetRangeUser(0.,0.5);

    TLegend *leg_shower_weight_energies_n;
    if(transformed == true) leg_shower_weight_energies_n= new TLegend(0.12,0.60,0.32,0.85);
    else leg_shower_weight_energies_n= new TLegend(0.12,0.60,0.32,0.85);
    leg_shower_weight_energies_n->SetTextSize(0.035);
    leg_shower_weight_energies_n->SetTextFont(42);
    leg_shower_weight_energies_n->AddEntry(g_shower_weight_n_2,"2 GeV","ap");
    leg_shower_weight_energies_n->AddEntry(g_shower_weight_n_4,"4 GeV","ap");
    leg_shower_weight_energies_n->AddEntry(g_shower_weight_n_6,"6 GeV","ap");
    leg_shower_weight_energies_n->AddEntry(g_shower_weight_n_10,"10 GeV","ap");
    leg_shower_weight_energies_n->AddEntry(g_shower_weight_n_20,"20 GeV","ap");
    leg_shower_weight_energies_n->AddEntry(g_shower_weight_n_40,"40 GeV","ap");
    leg_shower_weight_energies_n->AddEntry(g_shower_weight_n_80,"80 GeV","ap");
    leg_shower_weight_energies_n->AddEntry(g_shower_weight_n_150,"150 GeV","ap");
    leg_shower_weight_energies_n->SetFillColor(0);
    leg_shower_weight_energies_n->SetLineColor(0);
    leg_shower_weight_energies_n->SetShadowColor(0);
    leg_shower_weight_energies_n->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_weight_energies_n->SaveAs("shower_weight_energies_n_"+savetrans+"_w_errors.eps");
      c_shower_weight_energies_n->SaveAs("shower_weight_energies_n_"+savetrans+"_w_errors.png");
    }

    // Shower max, start and end
    // Shower weight max
    TGraphErrors *g_shower_weight_max = new TGraphErrors(*energyaxis, *shower_weight_max, *zeros, *zeros);
    graph_setup_add(g_shower_weight_max, "Shower maximum (w. energy)", kBlue);

    TMultiGraph *mg_shower_weight_max = new TMultiGraph();
    mg_shower_weight_max->Add(g_shower_weight_max);

    auto c_shower_weight_max = new TCanvas("c_shower_weight_max", "c_shower_weight_max", 800, 800);
    mg_shower_weight_max->Draw("AP");
    mg_shower_weight_max->SetTitle("Shower max (w. energy) ("+partstring+")");

    if (transformed) mg_shower_weight_max->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_weight_max->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_weight_max->GetYaxis()->SetTitle("Layer");
    mg_shower_weight_max->GetYaxis()->SetTitleOffset(1.4);
    mg_shower_weight_max->GetYaxis()->SetRangeUser(0,50);

    TLegend *leg_shower_weight_max;
    if(transformed == true) leg_shower_weight_max= new TLegend(0.45,0.65,0.65,0.85);
    else leg_shower_weight_max= new TLegend(0.45,0.65,0.65,0.85);
    leg_shower_weight_max->SetTextSize(0.035);
    leg_shower_weight_max->SetTextFont(42);
    leg_shower_weight_max->AddEntry(g_shower_weight_max,"Shower max. (w. energy)","ap");
    leg_shower_weight_max->SetFillColor(0);
    leg_shower_weight_max->SetLineColor(0);
    leg_shower_weight_max->SetShadowColor(0);
    leg_shower_weight_max->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_weight_max->SaveAs("shower_weight_max_"+savetrans+".eps");
      c_shower_weight_max->SaveAs("shower_weight_max_"+savetrans+".png");
    }

    // Shower weight average
    TGraphErrors *g_shower_weight_average = new TGraphErrors(*energyaxis, *shower_weight_average, *zeros, *zeros);
    graph_setup_add(g_shower_weight_average, "Shower average (w. energy)", kBlue);

    TMultiGraph *mg_shower_weight_average = new TMultiGraph();
    mg_shower_weight_average->Add(g_shower_weight_average);

    auto c_shower_weight_average = new TCanvas("c_shower_weight_average", "c_shower_weight_average", 800, 800);
    mg_shower_weight_average->Draw("AP");
    mg_shower_weight_average->SetTitle("Shower average (w. energy) ("+partstring+")");

    if (transformed) mg_shower_weight_average->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_weight_average->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_weight_average->GetYaxis()->SetTitle("Layer");
    mg_shower_weight_average->GetYaxis()->SetTitleOffset(1.4);
    mg_shower_weight_average->GetYaxis()->SetRangeUser(0,50);

    TLegend *leg_shower_weight_average;
    if(transformed == true) leg_shower_weight_average= new TLegend(0.45,0.65,0.65,0.85);
    else leg_shower_weight_average= new TLegend(0.45,0.65,0.65,0.85);
    leg_shower_weight_average->SetTextSize(0.035);
    leg_shower_weight_average->SetTextFont(42);
    leg_shower_weight_average->AddEntry(g_shower_weight_average,"Shower average (w. energy)","ap");
    leg_shower_weight_average->SetFillColor(0);
    leg_shower_weight_average->SetLineColor(0);
    leg_shower_weight_average->SetShadowColor(0);
    leg_shower_weight_average->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_weight_average->SaveAs("shower_weight_average_"+savetrans+".eps");
      c_shower_weight_average->SaveAs("shower_weight_average_"+savetrans+".png");
    }

    // Shower weight max layer
    TGraphErrors *g_shower_weight_max_layer = new TGraphErrors(*energyaxis, *shower_weight_max_layer, *zeros, *zeros);
    graph_setup_add(g_shower_weight_max_layer, "Shower maximum (w. energy)", kBlue);

    TMultiGraph *mg_shower_weight_max_layer = new TMultiGraph();
    mg_shower_weight_max_layer->Add(g_shower_weight_max_layer);

    auto c_shower_weight_max_layer = new TCanvas("c_shower_weight_max_layer", "c_shower_weight_max_layer", 800, 800);
    mg_shower_weight_max_layer->Draw("AP");
    mg_shower_weight_max_layer->SetTitle("W. energy shower max layer (avg.) ("+partstring+")");

    if (transformed) mg_shower_weight_max_layer->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_weight_max_layer->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_weight_max_layer->GetYaxis()->SetTitle("Layer");
    mg_shower_weight_max_layer->GetYaxis()->SetTitleOffset(1.4);
    mg_shower_weight_max_layer->GetYaxis()->SetRangeUser(-1,15);

    TLegend *leg_shower_weight_max_layer;
    if(transformed == true) leg_shower_weight_max_layer= new TLegend(0.45,0.65,0.65,0.85);
    else leg_shower_weight_max_layer= new TLegend(0.45,0.65,0.65,0.85);
    leg_shower_weight_max_layer->SetTextSize(0.035);
    leg_shower_weight_max_layer->SetTextFont(42);
    leg_shower_weight_max_layer->AddEntry(g_shower_weight_max_layer,"Shower max layer (weighted energy)","ap");
    leg_shower_weight_max_layer->SetFillColor(0);
    leg_shower_weight_max_layer->SetLineColor(0);
    leg_shower_weight_max_layer->SetShadowColor(0);
    leg_shower_weight_max_layer->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_weight_max_layer->SaveAs("shower_weight_max_layer_"+savetrans+".eps");
      c_shower_weight_max_layer->SaveAs("shower_weight_max_layer_"+savetrans+".png");
    }

    // Shower weight start layer
    TGraphErrors *g_shower_weight_start_layer = new TGraphErrors(*energyaxis, *shower_weight_start_layer, *zeros, *zeros);
    graph_setup_add(g_shower_weight_start_layer, "Shower start (w. energy)", kBlue);

    TMultiGraph *mg_shower_weight_start_layer = new TMultiGraph();
    mg_shower_weight_start_layer->Add(g_shower_weight_start_layer);

    auto c_shower_weight_start_layer = new TCanvas("c_shower_weight_start_layer", "c_shower_weight_start_layer", 800, 800);
    mg_shower_weight_start_layer->Draw("AP");
    mg_shower_weight_start_layer->SetTitle("W. energy shower start layer (avg.) ("+partstring+")");

    if (transformed) mg_shower_weight_start_layer->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_weight_start_layer->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_weight_start_layer->GetYaxis()->SetTitle("Layer");
    mg_shower_weight_start_layer->GetYaxis()->SetTitleOffset(1.4);
    mg_shower_weight_start_layer->GetYaxis()->SetRangeUser(-1,15);

    TLegend *leg_shower_weight_start_layer;
    if(transformed == true) leg_shower_weight_start_layer= new TLegend(0.45,0.65,0.65,0.85);
    else leg_shower_weight_start_layer= new TLegend(0.45,0.65,0.65,0.85);
    leg_shower_weight_start_layer->SetTextSize(0.035);
    leg_shower_weight_start_layer->SetTextFont(42);
    leg_shower_weight_start_layer->AddEntry(g_shower_weight_start_layer,"Shower start layer (w. energy)","ap");
    leg_shower_weight_start_layer->SetFillColor(0);
    leg_shower_weight_start_layer->SetLineColor(0);
    leg_shower_weight_start_layer->SetShadowColor(0);
    leg_shower_weight_start_layer->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_weight_start_layer->SaveAs("shower_weight_start_layer_"+savetrans+".eps");
      c_shower_weight_start_layer->SaveAs("shower_weight_start_layer_"+savetrans+".png");
    }


    // Shower weight start end
    TGraphErrors *g_shower_weight_end_layer = new TGraphErrors(*energyaxis, *shower_weight_end_layer, *zeros, *zeros);
    graph_setup_add(g_shower_weight_end_layer, "Shower end (w. energy)", kBlue);

    TMultiGraph *mg_shower_weight_end_layer = new TMultiGraph();
    mg_shower_weight_end_layer->Add(g_shower_weight_end_layer);

    auto c_shower_weight_end_layer = new TCanvas("c_shower_weight_end_layer", "c_shower_weight_end_layer", 800, 800);
    mg_shower_weight_end_layer->Draw("AP");
    mg_shower_weight_end_layer->SetTitle("W. energy shower end layer (avg.) ("+partstring+")");

    if (transformed) mg_shower_weight_end_layer->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_shower_weight_end_layer->GetXaxis()->SetTitle("E[GeV]");
    mg_shower_weight_end_layer->GetYaxis()->SetTitle("Layer");
    mg_shower_weight_end_layer->GetYaxis()->SetTitleOffset(1.4);
    mg_shower_weight_end_layer->GetYaxis()->SetRangeUser(-1,15);

    TLegend *leg_shower_weight_end_layer;
    if(transformed == true) leg_shower_weight_end_layer= new TLegend(0.45,0.65,0.65,0.85);
    else leg_shower_weight_end_layer= new TLegend(0.45,0.65,0.65,0.85);
    leg_shower_weight_end_layer->SetTextSize(0.035);
    leg_shower_weight_end_layer->SetTextFont(42);
    leg_shower_weight_end_layer->AddEntry(g_shower_weight_end_layer,"Shower end layer (w. energy)","ap");
    leg_shower_weight_end_layer->SetFillColor(0);
    leg_shower_weight_end_layer->SetLineColor(0);
    leg_shower_weight_end_layer->SetShadowColor(0);
    leg_shower_weight_end_layer->Draw();
    addCaliceLogo();
    if(save == true){
      c_shower_weight_end_layer->SaveAs("shower_weight_end_layer_"+savetrans+".eps");
      c_shower_weight_end_layer->SaveAs("shower_weight_end_layer_"+savetrans+".png");
    }


}
