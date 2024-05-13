#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include "TGraph.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TLeaf.h"
#include "TH1F.h"
#include "TCanvas.h"

#define N_BINS 40

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

void addProduction(){
  TPad *p1 = new TPad("img", "img", 0.1, 0.905, 0.5, 0.925);
  p1->Draw();
  p1->cd();
  TString title = "PID_100_GeV_e_pi_mu";
  TText* t = new TText(0.01,0.1,title);
  t->SetTextColor(kRed);
  t->SetTextSize(0.99);
  t->Draw();
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

void findTH1values(TH1F * hA, TH1F * hB, TH1F * hC, float &x_min, float &x_max, float &y_max) {
  float max_y_A = hA->GetBinContent(hA->GetMaximumBin());
  float max_y_B = hB->GetBinContent(hB->GetMaximumBin());
  float max_y_C = hC->GetBinContent(hC->GetMaximumBin());

  //cout<<"Max (e, mu, pi): ("<<max_y_A<<", "<<max_y_B<<", "<<max_y_C<<")"<<endl;

  int b_min_A = hA->FindFirstBinAbove(0,1);
  int b_min_B = hB->FindFirstBinAbove(0,1);
  int b_min_C = hC->FindFirstBinAbove(0,1);

  int b_max_A = hA->FindLastBinAbove(0,1);
  int b_max_B = hB->FindLastBinAbove(0,1);
  int b_max_C = hC->FindLastBinAbove(0,1);

  int b_min = min(b_min_A, min(b_min_B,b_min_C));
  int b_max = max(b_max_A, max(b_max_B,b_max_C));

  x_min = hA->GetBinCenter(b_min-1);
  x_max = hA->GetBinCenter(b_max+1);

  y_max = 1.01*max(max_y_A, max(max_y_B,max_y_C));
  //cout<<"(min x, max x, max y) = ("<<x_min<<", "<<x_max<<", "<<y_max<<")"<<endl;

  return 0;
}

void drawLikelihood(TString title, TH1F * hA, TH1F * hB, TH1F * hC, bool save = false) {
  float x_min = 0.;
  float x_max = 0.;
  float y_max = 0.;
 
  findTH1values(hA, hB, hC, x_min, x_max, y_max);
  hA->GetYaxis()->SetMaxDigits(3);
  hA->SetTitle(title);
  hA->SetXTitle(title);
  hA->SetYTitle("entries (%)");

  float sum_A = 0.;
  float sum_B = 0.;
  float sum_C = 0.;
  for(int i=1; i<N_BINS+1; i++) {
    sum_A += hA->GetBinContent(i);
    sum_B += hB->GetBinContent(i);
    sum_C += hC->GetBinContent(i);
  }
  //cout<<"Sums: "<<sum_A<<" "<<sum_B<<" "<<sum_C<<endl;
  //cout<<"Integrals: "<<hA->Integral()<<" "<<hB->Integral()<<" "<<hC->Integral()<<endl;
  
  auto c = new TCanvas("c_PID_100_GeV_e_pi_mu_"+title, "c_PID_100_GeV_e_pi_mu_"+title, 800, 800);
  c->cd();
  hA->SetLineColor(kGray+2);
  hB->SetLineColor(kCyan-3);
  hC->SetLineColor(kRed-4);
  hA->SetLineWidth(3);
  hB->SetLineWidth(3);
  hC->SetLineWidth(3);
  hA->Scale(100./sum_A);
  hB->Scale(100./sum_B);
  hC->Scale(100./sum_C);
  //cout<<"y max: "<<y_max<<" %: "<<100.*y_max/sum_A<<endl;
  hA->GetYaxis()->SetRangeUser(0,100.*y_max/sum_A);
  hA->GetYaxis()->SetTitleOffset(1.4);
  hA->Draw("histo");
  hB->Draw("histosame");
  hC->Draw("histosame");
  
  float xleg=0.5;
  if(title=="e-likelihood") xleg = 0.3;
  if(title=="#pi-likelihood") xleg = 0.55;
  if(title=="#mu-likelihood") xleg = 0.55;
  TLegend *leg;
  leg= new TLegend(xleg,0.70,xleg+0.2,0.85);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->AddEntry(hA,"e-","l");
  leg->AddEntry(hB,"#pi-","l");
  leg->AddEntry(hC,"#mu-","l");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->Draw();
  addProduction();
  c->cd();
  addCaliceLogo(true);

  if(save == true){
    c->SaveAs("c_PID_100_GeV_e_pi_mu_"+title+".eps");
    c->SaveAs("c_PID_100_GeV_e_pi_mu_"+title+".png");
  }

  // Building the efficiency histos
  double cut_axis[N_BINS+1];
  for(int i=0; i<N_BINS+1; i++){
    cut_axis[i] = 0. + i*0.0250;
    //cout<<cut_axis[i]<<endl;
  }
  double hA_c[N_BINS+1], hB_c[N_BINS+1], hC_c[N_BINS+1];
  double h_001_x[1], h_01_x[1], h_1_x[1], h_10_x[1];
  double h_001_y[1], h_01_y[1], h_1_y[1], h_10_y[1];
  h_001_x[0] = -1.;
  h_01_x[0] = -1.;
  h_1_x[0] = -1.;
  h_10_x[0] = -1.;
  h_001_y[0] = -1.;
  h_01_y[0] = -1.;
  h_1_y[0] = -1.;
  h_10_y[0] = -1.;
  bool found_001 = false;
  bool found_01 = false;
  bool found_1 = false;
  bool found_10 = false;

  float integral_A = 0.;
  float integral_B = 0.;
  float integral_C = 0.;

  for(int i=0;i<N_BINS+1;i++){
    if(i==0){
      hA_c[i] = 100.;
      hB_c[i] = 100.;
      hC_c[i] = 100.;
      continue;
    }
    float A_content = hA->GetBinContent(i);
    float B_content = hB->GetBinContent(i);
    float C_content = hC->GetBinContent(i);
    integral_A += A_content;
    integral_B += B_content;
    integral_C += C_content;
    //cout<<title<<" A (bin "<<i+1<<") "<<A_content<<" cumulated: "<<100.-integral_A<<endl;
    //cout<<title<<" B (bin "<<i+1<<") "<<B_content<<" cumulated: "<<100.-integral_B<<endl;
    //cout<<title<<" C (bin "<<i+1<<") "<<C_content<<" cumulated: "<<100.-integral_C<<endl;
    hA_c[i] = 100.-integral_A;
    hB_c[i] = 100.-integral_B;
    hC_c[i] = 100.-integral_C;
    if(hB_c[i] < 0.) hB_c[i] = 0.;
    if(hC_c[i] < 0.) hC_c[i] = 0.;
  }
      
  auto hA_cumulative = new TGraph(N_BINS+1,cut_axis,hA_c);
  auto hB_cumulative = new TGraph(N_BINS+1,cut_axis,hB_c);
  auto hC_cumulative = new TGraph(N_BINS+1,cut_axis,hC_c);
  
  int NPoints = 1000;
  float tol = 0.01;
  for(int i=1; i<NPoints; i++) {
    float x = float(i)/float(NPoints);
    float A_content = hA_cumulative->Eval(x);
    float B_content = hB_cumulative->Eval(x);
    float C_content = hC_cumulative->Eval(x);
    //cout<<"x: "<<x<<". (A, B, C) = ("<<A_content<<", "<<B_content<<", "<<C_content<<")"<<endl;
    float A_impurity = 100.*(B_content + C_content)/(A_content + B_content + C_content);
    float B_impurity = 100.*(A_content + C_content)/(A_content + B_content + C_content);
    float C_impurity = 100.*(A_content + B_content)/(A_content + B_content + C_content);
    //cout<<"x: "<<x<<". Impurity (A, B, C) = "<<A_impurity<<", "<<B_impurity<<", "<<C_impurity<<")"<<endl;
    if( title=="e-likelihood" ){
      if( A_impurity < 0.01*(1+tol) ){
	h_001_y[0] = A_content;
	h_001_x[0] = x;
	found_001 = true;
	break;
      }
    }
    if( title=="#pi-likelihood" ){
      if( B_impurity < 0.01*(1+tol) ){
        h_001_y[0] = B_content;
        h_001_x[0] = x;
        found_001 = true;
        break;
      }
    }
    if( title=="#mu-likelihood" ){
      if( C_impurity < 0.01*(1+tol) ){
        h_001_y[0] = C_content;
        h_001_x[0] = x;
        found_001 = true;
        break;
      }
    }
  }
  for(int i=1; i<NPoints; i++) {
    float x = float(i)/float(NPoints);
    float A_content = hA_cumulative->Eval(x);
    float B_content = hB_cumulative->Eval(x);
    float C_content = hC_cumulative->Eval(x);
    //cout<<"x: "<<x<<". (A, B, C) = ("<<A_content<<", "<<B_content<<", "<<C_content<<")"<<endl;
    float A_impurity = 100.*(B_content + C_content)/(A_content + B_content + C_content);
    float B_impurity = 100.*(A_content + C_content)/(A_content + B_content + C_content);
    float C_impurity = 100.*(A_content + B_content)/(A_content + B_content + C_content);
    //cout<<"impurity (A, B, C) = "<<A_impurity<<", "<<B_impurity<<", "<<C_impurity<<")"<<endl;
    if( title=="e-likelihood" ){
      if( A_impurity < 0.1*(1+tol) ){
        h_01_y[0] = A_content;
        h_01_x[0] = x;
        found_01 = true;
        break;
      }
    }
    if( title=="#pi-likelihood" ){
      if( B_impurity < 0.1*(1+tol) ){
        h_01_y[0] = B_content;
        h_01_x[0] = x;
        found_01 = true;
        break;
      }
    }
    if( title=="#mu-likelihood" ){
      if( C_impurity < 0.1*(1+tol) ){
        h_01_y[0] = C_content;
        h_01_x[0] = x;
        found_01 = true;
        break;
      }
    }
  }
  for(int i=1; i<NPoints; i++) {
    float x = float(i)/float(NPoints);
    float A_content = hA_cumulative->Eval(x);
    float B_content = hB_cumulative->Eval(x);
    float C_content = hC_cumulative->Eval(x);
    //cout<<"x: "<<x<<". (A, B, C) = ("<<A_content<<", "<<B_content<<", "<<C_content<<")"<<endl;
    float A_impurity = 100.*(B_content + C_content)/(A_content + B_content + C_content);
    float B_impurity = 100.*(A_content + C_content)/(A_content + B_content + C_content);
    float C_impurity = 100.*(A_content + B_content)/(A_content + B_content + C_content);
    //cout<<"impurity (A, B, C) = "<<A_impurity<<", "<<B_impurity<<", "<<C_impurity<<")"<<endl;
    if( title=="e-likelihood" ){
      if( A_impurity < 1.*(1+tol) ){
        h_1_y[0] = A_content;
        h_1_x[0] = x;
        found_1 = true;
        break;
      }
    }
    if( title=="#pi-likelihood" ){
      if( B_impurity < 1.*(1+tol) ){
        h_1_y[0] = B_content;
        h_1_x[0] = x;
        found_1 = true;
        break;
      }
    }
    if( title=="#mu-likelihood" ){
      if( C_impurity < 1.*(1+tol) ){
        h_1_y[0] = C_content;
        h_1_x[0] = x;
        found_1 = true;
        break;
      }
    }
  }
  for(int i=1; i<NPoints; i++) {
    float x = float(i)/float(NPoints);
    float A_content = hA_cumulative->Eval(x);
    float B_content = hB_cumulative->Eval(x);
    float C_content = hC_cumulative->Eval(x);
    //cout<<"x: "<<x<<". (A, B, C) = ("<<A_content<<", "<<B_content<<", "<<C_content<<")"<<endl;
    float A_impurity = 100.*(B_content + C_content)/(A_content + B_content + C_content);
    float B_impurity = 100.*(A_content + C_content)/(A_content + B_content + C_content);
    float C_impurity = 100.*(A_content + B_content)/(A_content + B_content + C_content);
    //cout<<"impurity (A, B, C) = "<<A_impurity<<", "<<B_impurity<<", "<<C_impurity<<")"<<endl;
    if( title=="e-likelihood" ){
      if( A_impurity < 10.*(1+tol) ){
        h_10_y[0] = A_content;
        h_10_x[0] = x;
        found_10 = true;
        break;
      }
    }
    if( title=="#pi-likelihood" ){
      if( B_impurity < 10.*(1+tol) ){
        h_10_y[0] = B_content;
        h_10_x[0] = x;
        found_10 = true;
        break;
      }
    }
    if( title=="#mu-likelihood" ){
      if( C_impurity < 10.*(1+tol) ){
        h_10_y[0] = C_content;
        h_10_x[0] = x;
        found_10 = true;
        break;
      }
    }
  }


  auto h_001 = new TGraph(1,h_001_x,h_001_y);
  auto h_01 = new TGraph(1,h_01_x,h_01_y);
  auto h_1 = new TGraph(1,h_1_x,h_1_y);
  auto h_10 = new TGraph(1,h_10_x,h_10_y);
  
  for(int i=0;i<N_BINS+1;i++) cout<<hA_c[i]<<" "<<hB_c[i]<<" "<<hC_c[i]<<endl;
  cout<<"P(0.1) = ("<<h_01_x[0]<<","<<h_01_y[0]<<"); "<<"P(1) = ("<<h_1_x[0]<<","<<h_1_y[0]<<"); "<<"P(10) = ("<<h_10_x[0]<<","<<h_10_y[0]<<");"<<endl;
  
  hA_cumulative->GetYaxis()->SetRangeUser(0,105.);
  hA_cumulative->GetXaxis()->SetRangeUser(-0.05,1.05);
  hA_cumulative->GetYaxis()->SetMaxDigits(3);
  hA_cumulative->SetTitle("Efficiencies for "+title+" cut;"+title+" cut;efficiency (%)");
  hA_cumulative->GetYaxis()->SetTitleOffset(1.4);

  auto c2 = new TCanvas("c_PID_100_GeV_e_pi_mu_"+title+"_cut", "c_PID_100_GeV_e_pi_mu_"+title+"_cut", 800, 800);
  c2->cd();
  hA_cumulative->SetLineColor(kGray+2);
  hB_cumulative->SetLineColor(kRed-4);
  hC_cumulative->SetLineColor(kCyan-3);
  hA_cumulative->SetLineWidth(3);
  hB_cumulative->SetLineWidth(3);
  hC_cumulative->SetLineWidth(3);
  hA_cumulative->Draw("ALP");
  hB_cumulative->Draw("LP");
  hC_cumulative->Draw("LP");
  h_001->SetMarkerStyle(33);
  h_001->SetMarkerSize(2);
  h_001->SetMarkerColor(77);
  h_01->SetMarkerStyle(22);
  h_01->SetMarkerSize(2);
  h_01->SetMarkerColor(92);
  h_1->SetMarkerStyle(20);
  h_1->SetMarkerSize(2);
  h_1->SetMarkerColor(94);
  h_10->SetMarkerStyle(21);
  h_10->SetMarkerSize(2);
  h_10->SetMarkerColor(98);

  if(found_10 == true) h_10->Draw("P");
  if(found_1 == true) h_1->Draw("P");
  if(found_01 == true) h_01->Draw("P");
  if(found_001 == true) h_001->Draw("P");

  float xleg_c = 0.6;
  float yleg_c = 0.6;
  if(title=="e-likelihood") xleg_c = 0.3;
  if(title=="#pi-likelihood") {xleg_c = 0.65; yleg_c = 0.7;}
  if(title=="#mu-likelihood") {xleg_c = 0.65; yleg_c = 0.7;}
  TLegend *leg_c;
  leg_c= new TLegend(xleg_c,yleg_c,xleg_c+0.2,yleg_c+0.15);
  leg_c->SetTextSize(0.025);
  leg_c->SetTextFont(42);
  leg_c->AddEntry(hA,"e-","l");
  leg_c->AddEntry(hB,"#pi-","l");
  leg_c->AddEntry(hC,"#mu-","l");
  leg_c->SetFillColor(0);
  leg_c->SetLineColor(0);
  leg_c->SetShadowColor(0);
  leg_c->Draw();
  
  if( (found_001 == true) or (found_01 == true) or (found_1 == true) or (found_10 == true)) {
    float xleg_WP = 0.6;
    float yleg_WP = 0.4;
    if(title=="e-likelihood") xleg_WP = 0.3;
    if(title=="#pi-likelihood") {xleg_WP = 0.65; yleg_WP = 0.5;}
    if(title=="#mu-likelihood") {xleg_WP = 0.65; yleg_WP = 0.5;}
    TLegend *leg_WP;
    leg_WP= new TLegend(xleg_WP,yleg_WP,xleg_WP+0.2,yleg_WP+0.15);
    leg_WP->SetTextSize(0.025);
    leg_WP->SetTextFont(42);
    if(found_001 == true) leg_WP->AddEntry(h_001,"99.99% Purity","P");
    if(found_01 == true) leg_WP->AddEntry(h_01,"99.9% Purity","P");
    if(found_1 == true) leg_WP->AddEntry(h_1,"99% Purity","P");
    if(found_10 == true) leg_WP->AddEntry(h_10,"90% Purity","P");
    leg_WP->SetFillColor(0);
    leg_WP->SetLineColor(0);
    leg_WP->SetShadowColor(0);
    leg_WP->Draw();
  }
  addProduction();
  c2->cd();
  addCaliceLogo(true);

  if(save == true){
    c2->SaveAs("c_PID_100_GeV_e_pi_mu_"+title+"_cut.eps");
    c2->SaveAs("c_PID_100_GeV_e_pi_mu_"+title+"_cut.png");
  }

}

void drawVariablesLegend(vector<TString> varnames) {
  TLatex text_vars;
  text_vars.SetTextSize(0.27);
  float linesize = 0.;
  float firstline = 0.;
  float secondline = 0.;

  for(int i=0;i<varnames.size();i++){
    //cout<<"string size: "<<varnames.at(i).Length()<<", ";
    if(i==0) {
      text_vars.DrawLatex(0.001,0.8,"V"+TString::Format("%i:",i+1)+varnames.at(i));
    }
    else{
      if(linesize < 80) {
        text_vars.DrawLatex(0.001+linesize*0.011,0.8,"V"+TString::Format("%i:",i+1)+varnames.at(i));
	firstline = linesize + float(varnames.at(i).Length())+3.+2.;
      }
      else if( (linesize > 79) and (linesize < 160) ) {
        text_vars.DrawLatex(0.001+(linesize-firstline)*0.011,0.5,"V"+TString::Format("%i:",i+1)+varnames.at(i));
        secondline = linesize + float(varnames.at(i).Length())+3.+2.;
      }
      else if(linesize > 159) {
        text_vars.DrawLatex(0.001+(linesize-secondline)*0.011,0.2,"V"+TString::Format("%i:",i+1)+varnames.at(i));
      }
    }
    linesize += float(varnames.at(i).Length())+3.+2.;
    //cout<<"linesize: "<<linesize<<endl;
  }
}

void drawCorrelation(TString particle, vector<TString> varnames, TH2F * hA, bool save=false) {

  auto c_A = new TCanvas("c_PID_100_GeV_e_pi_mu_X_Correlation_"+particle,"c_PID_100_GeV_e_pi_mu_X_Correlation_"+particle,800,800);
  c_A->cd();
  hA->SetTitle("|Variable correlation| ("+particle+")");
  hA->GetXaxis()->LabelsOption("h");
  hA->GetXaxis()->SetLabelSize(0.04);
  hA->GetXaxis()->SetLabelOffset(0.005);
  hA->GetXaxis()->SetTickSize(0.);
  hA->GetYaxis()->SetLabelSize(0.04);
  hA->GetYaxis()->SetLabelOffset(0.005);
  hA->GetYaxis()->SetTickSize(0.);
  hA->GetZaxis()->SetLabelSize(0.025);
  hA->SetTitleOffset(0.8);
  hA->Draw("colz");
  addCaliceLogo(true);
  c_A->cd();
  TPad *pad = new TPad("padLegend","",0.1,0,0.9,0.06,5);
  pad->SetFillStyle(4000);
  pad->Draw();
  pad->cd();
  drawVariablesLegend(varnames);
  c_A->cd();
  addProduction();
  
  if(save == true){
    c_A->SaveAs("c_PID_100_GeV_e_pi_mu_correlation_"+particle+".eps");
    c_A->SaveAs("c_PID_100_GeV_e_pi_mu_correlation_"+particle+".png");
  }
}

void drawCorrelationPlots(TH2F * hA, TH2F * hB, TH2F * hC, bool save = false) {

  Int_t nbins_x = hA->GetNbinsX();
  Int_t nbins_y = hA->GetNbinsY();
  vector<TString> varnames;
  for(int i = 1; i<nbins_x+1; i++) varnames.push_back(hA->GetXaxis()->GetLabels()->At(i-1)->GetName());

  TH2F * A_cor = new TH2F("e_correlation","e sample correlation between variables",nbins_x,-0.5,nbins_x+0.5,nbins_y,-0.5,nbins_y+0.5);
  TH2F * B_cor = new TH2F("#pi_correlation","#pi sample correlation between variables",nbins_x,-0.5,nbins_x+0.5,nbins_y,-0.5,nbins_y+0.5);
  TH2F * C_cor = new TH2F("#mu_correlation","#mu sample correlation between variables",nbins_x,-0.5,nbins_x+0.5,nbins_y,-0.5,nbins_y+0.5);
  TH2F * max_correlation = new TH2F("max_correlation","Maximum correlation between variables",nbins_x,-0.5,nbins_x+0.5,nbins_y,-0.5,nbins_y+0.5);
  for(int i = 1; i<nbins_x+1; i++){
    for(int j = 1; j<nbins_y+1; j++){
      float content_A = abs(hA->GetBinContent(i,j));
      float content_B = abs(hB->GetBinContent(i,j));
      float content_C = abs(hC->GetBinContent(i,j));
      if(content_A == 0.) content_A = 0.001;
      if(content_B == 0.) content_B = 0.001;
      if(content_C == 0.) content_C = 0.001;
      float maxvalue = max(content_A,max(content_B,content_C));
      max_correlation->SetBinContent(i,j,maxvalue);
      A_cor->SetBinContent(i,j,content_A);
      B_cor->SetBinContent(i,j,content_B);
      C_cor->SetBinContent(i,j,content_C);
    }
    A_cor->GetXaxis()->SetBinLabel(i,"V"+TString::Format("%i",i));
    A_cor->GetYaxis()->SetBinLabel(i,"V"+TString::Format("%i",i));
    B_cor->GetXaxis()->SetBinLabel(i,"V"+TString::Format("%i",i));
    B_cor->GetYaxis()->SetBinLabel(i,"V"+TString::Format("%i",i));
    C_cor->GetXaxis()->SetBinLabel(i,"V"+TString::Format("%i",i));
    C_cor->GetYaxis()->SetBinLabel(i,"V"+TString::Format("%i",i));
    max_correlation->GetXaxis()->SetBinLabel(i,"V"+TString::Format("%i",i));
    max_correlation->GetYaxis()->SetBinLabel(i,"V"+TString::Format("%i",i));
  }

  drawCorrelation("e-", varnames, A_cor, save);
  drawCorrelation("#pi-", varnames, B_cor, save);
  drawCorrelation("#mu-", varnames, C_cor, save);
  drawCorrelation("max", varnames, max_correlation, save);
}

void Performance_plots(bool save = false){

  TFile* outputfile=TFile::Open("../ROOTFiles/PID_100_GeV_e_pi_mu_TMVA.root");
  //outputfile.cd("dataloader");
  //TTree* TestTree=(TTree*) outputfile->Get("dataloader/TestTree");
  //Int_t test_entries=TestTree->GetEntries();
  //cout<<"Test TTree entries: "<<test_entries<<endl;
  
  // First we'll get the likelihood histos
  //TH1F* Likelihood_e_for_e=new TH1F("likelihood_e_for_e","likelihood_e_for_e",40,0,1);
  outputfile->cd();
  
  TH1F* A_likelihood_for_A = (TH1F*) outputfile->Get("dataloader/Method_BDT/myMVA/MVA_myMVA_Test_Sample_A_prob_for_Sample_A");
  TH1F* A_likelihood_for_B = (TH1F*) outputfile->Get("dataloader/Method_BDT/myMVA/MVA_myMVA_Test_Sample_B_prob_for_Sample_A");
  TH1F* A_likelihood_for_C = (TH1F*) outputfile->Get("dataloader/Method_BDT/myMVA/MVA_myMVA_Test_Sample_C_prob_for_Sample_A");
  
  TH1F* B_likelihood_for_A = (TH1F*) outputfile->Get("dataloader/Method_BDT/myMVA/MVA_myMVA_Test_Sample_A_prob_for_Sample_B");
  TH1F* B_likelihood_for_B = (TH1F*) outputfile->Get("dataloader/Method_BDT/myMVA/MVA_myMVA_Test_Sample_B_prob_for_Sample_B");
  TH1F* B_likelihood_for_C = (TH1F*) outputfile->Get("dataloader/Method_BDT/myMVA/MVA_myMVA_Test_Sample_C_prob_for_Sample_B");
  
  TH1F* C_likelihood_for_A = (TH1F*) outputfile->Get("dataloader/Method_BDT/myMVA/MVA_myMVA_Test_Sample_A_prob_for_Sample_C");
  TH1F* C_likelihood_for_B = (TH1F*) outputfile->Get("dataloader/Method_BDT/myMVA/MVA_myMVA_Test_Sample_B_prob_for_Sample_C");
  TH1F* C_likelihood_for_C = (TH1F*) outputfile->Get("dataloader/Method_BDT/myMVA/MVA_myMVA_Test_Sample_C_prob_for_Sample_C");


  TH2F* correlation_matrix_A = (TH2F*) outputfile->Get("dataloader/CorrelationMatrixSample_A");
  TH2F* correlation_matrix_B = (TH2F*) outputfile->Get("dataloader/CorrelationMatrixSample_B");
  TH2F* correlation_matrix_C = (TH2F*) outputfile->Get("dataloader/CorrelationMatrixSample_C");
  
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleStyle(0);
  //gStyle->SetTitleX(0.2);
  gStyle->SetMarkerSize(0.2);
  
  drawLikelihood("e-likelihood", A_likelihood_for_A, A_likelihood_for_B, A_likelihood_for_C, save);
  drawLikelihood("#pi-likelihood", B_likelihood_for_A, B_likelihood_for_B, B_likelihood_for_C, save);
  drawLikelihood("#mu-likelihood", C_likelihood_for_A, C_likelihood_for_B, C_likelihood_for_C, save);
 
  drawCorrelationPlots(correlation_matrix_A, correlation_matrix_B, correlation_matrix_C, save);
}
