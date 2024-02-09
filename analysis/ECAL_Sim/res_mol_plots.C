#define N_ENERGIES 18

void addCaliceLogo(bool WIP = true){
  TImage *img = TImage::Open("style/CALICELogo_18pc.png");
  img->SetConstRatio(kTRUE);
  img->SetImageCompression(0);
  TPad *p1 = new TPad("img", "img", 0.80, 0.90, 1.0, 1.0);
  p1->Draw();
  p1->cd();
  img->Draw();

  if(WIP == true){
    TPad *p2 = new TPad("img", "img", 0.025, 0.01, 1.0, 0.21);
    p1->cd();
    p2->Draw();
    p2->cd();
    TText* t = new TText(0.1,0.3,"Work in progress");
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

void res_mol_plots(string filename, bool transformed = true, bool save = false){
    
    TFile *file = new TFile(filename.c_str(), "read");
    // TVectorD * energies = file->GetObject("energies");
    TString savetrans;
    if(transformed == false) savetrans="Linear";
    else savetrans="InvSq";
    
    TVectorD *energies;
    TVectorD *energies_tr;
    
    TVectorD *res_nhit, *mu_nhit, *sig_nhit;
    TVectorD *res_nhit_masked, *mu_nhit_masked, *sig_nhit_masked;
    TVectorD *res_sume, *mu_sume, *sig_sume;
    TVectorD *res_sume_masked, *mu_sume_masked, *sig_sume_masked;
    TVectorD *res_weight, *mu_weight, *sig_weight;
    TVectorD *res_weight_masked, *mu_weight_masked, *sig_weight_masked;
    
    TVectorD *mol, *mol_masked, *mol_sig, *mol_sig_masked;

    TVectorD *zeros = new TVectorD(N_ENERGIES);
    
    double zerovalues[N_ENERGIES] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}; 
    TVectorD zerovectorref(N_ENERGIES, zerovalues);
    zeros = &zerovectorref;
    
    //zeros = new TVectorD(N_ENERGIES); 
    //    for (int j = 0; j < N_ENERGIES; j++) zeros[j] = 0.;
    
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
    
    file->GetObject("mol", mol);
    file->GetObject("mol_masked", mol_masked);
    file->GetObject("mol_sig", mol_sig);
    file->GetObject("mol_sig_masked", mol_sig_masked);

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
    
    
    /*
    // Moliere
    // TGraphErrors *g_mol = new TGraphErrors(N_ENERGIES, energies_tr, mol, zeros, mol_sig);   
    TGraphErrors *g_mol = new TGraphErrors(*energies_tr, *mol, zeros, *mol_sig);   
    graph_setup_add(g_mol, "Moliere Radius", kBlue);
    
    // TGraphErrors *g_mol_masked = new TGraphErrors(N_ENERGIES, energies_tr, mol_masked, zeros, mol_sig_masked);   
    TGraphErrors *g_mol_masked = new TGraphErrors(*energies_tr, *mol_masked, zeros, *mol_sig_masked);   
    graph_setup_add(g_mol_masked, "Moliere Radius (masked)", kRed);

    TMultiGraph *mg_mol = new TMultiGraph();
    mg_mol->Add(g_mol);
    // mg_mol->Add(g_mol_masked);

    auto c2 = new TCanvas("c2", "c2", 800, 600);
    
    mg_mol->Draw("AL");
    mg_mol->SetTitle("Moliere radius vs energy");
    // mg_mol->GetXaxis()->SetTitle("Energy [GeV]");
    if (transformed) mg_mol->GetXaxis()->SetTitle("1/#sqrt{E[GeV]}");
    else mg_mol->GetXaxis()->SetTitle("E[GeV]");
    mg_mol->GetYaxis()->SetTitle("Radius (95%) [mm]");
    if (!transformed){
        mg_mol->GetXaxis()->SetNdivisions(-502);
        int bin_index;
        for (int j = 0; j < N_ENERGIES; j++) {
            bin_index = mg_mol->GetXaxis()->FindBin(((*energies_tr))[j]);
            if (j == 8 || j == 7 || j == 9 ) continue;
            mg_mol->GetXaxis()->SetBinLabel(bin_index, to_string((int)round(((*energies))[j])).c_str());
            }
        mg_mol->GetXaxis()->LabelsOption("h");
        mg_mol->GetXaxis()->SetLabelSize(0.05);
    }
    
    c2->BuildLegend();
    c2->SaveAs("moliere_ecal_sim.eps");
    c2->SaveAs("moliere_ecal_sim.root");
    */
}
