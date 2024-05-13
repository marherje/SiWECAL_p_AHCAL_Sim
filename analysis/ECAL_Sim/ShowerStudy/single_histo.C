#define N_ECAL_LAYERS 15

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

void findTH1values(TH1 * he, TH1 * hmu, TH1 * hpi, float &x_min, float &x_max, float &y_max, bool &rescalable) {
  float max_x_e = he->GetBinContent(he->GetMaximumBin());
  float max_x_mu = hmu->GetBinContent(hmu->GetMaximumBin());
  float max_x_pi = hpi->GetBinContent(hpi->GetMaximumBin());
  
  cout<<"Max (e, mu, pi): ("<<max_x_e<<", "<<max_x_mu<<", "<<max_x_pi<<")"<<endl;

  //int entries = he->GetNbinsX();                                                                                                                                                                                                                                                            
  int b_min_e = he->FindFirstBinAbove(0,1);
  int b_min_mu = hmu->FindFirstBinAbove(0,1);
  int b_min_pi = hpi->FindFirstBinAbove(0,1);

  int b_max_e = he->FindLastBinAbove(0,1);
  int b_max_mu = hmu->FindLastBinAbove(0,1);
  int b_max_pi = hpi->FindLastBinAbove(0,1);

  int b_min = min(b_min_e, min(b_min_mu,b_min_pi));
  int b_max = max(b_max_e, max(b_max_mu,b_max_pi));

  if(b_max != N_ECAL_LAYERS + 1) rescalable = true;

  x_min = he->GetBinCenter(b_min-1);
  x_max = he->GetBinCenter(b_max+1);

  y_max = 1.01*max(max_x_e, max(max_x_mu,max_x_pi));
  cout<<"(min x, max x, max y) = ("<<x_min<<", "<<x_max<<", "<<y_max<<")"<<endl;

  return 0;
}

void drawHistosTH1(TString varname, TString energy, TH1 * he, TH1 * hmu, TH1 * hpi, bool save = false) {
  float x_min = 0.;
  float x_max = 0.;
  float y_max = 0.;
  bool rescalable = false;
  findTH1values(he, hmu, hpi, x_min, x_max, y_max, rescalable);
  he->GetYaxis()->SetRangeUser(0,y_max);
  he->GetYaxis()->SetMaxDigits(3);
  if(rescalable == true) he->GetXaxis()->SetRangeUser(x_min,x_max);
  /*
  TString power = "";
  if(y_max > 10000) {
    he->Scale(1./1000);
    hmu->Scale(1./1000);
    hpi->Scale(1./1000);
    he->GetYaxis()->SetRangeUser(0,y_max/1000);
    power = " (10^{3})";
  }
  else if(y_max > 1000) {
    he->Scale(1./100);
    hmu->Scale(1./100);
    hpi->Scale(1./100);
    he->GetYaxis()->SetRangeUser(0,y_max/100);
    power = " (10^{2})";
  }
  else if(y_max > 100) {
    he->Scale(1./10);
    hmu->Scale(1./10);
    hpi->Scale(1./10);
    he->GetYaxis()->SetRangeUser(0,y_max/10);
    power = " (10^{1})";
  }
  */
  //he->SetTitle(varname);
  he->SetXTitle("Value");
  he->SetYTitle("Entries");

  auto c = new TCanvas("c_"+varname+"_"+energy+"GeV", "c_"+varname+"_"+energy+"GeV", 800, 800);
  c->cd();
  if(y_max > 1000){
    he->SetMinimum(1);
    hmu->SetMinimum(1);
    hpi->SetMinimum(1);
    c->SetLogy();
  }
  he->SetLineColor(kGray+2);
  hmu->SetLineColor(kCyan-3);
  hpi->SetLineColor(kRed-4);
  he->SetLineWidth(3);
  hmu->SetLineWidth(3);
  hpi->SetLineWidth(3);

  he->Draw("histo");
  hmu->Draw("histosame");
  hpi->Draw("histosame");

  float xleg=0.;
  if( he->GetBinCenter(he->GetMaximumBin()) < (x_max-x_min)/2) xleg = 0.7;
  if( he->GetBinCenter(he->GetMaximumBin()) > (x_max-x_min)/2) xleg = 0.2;
  cout<<he->GetBinCenter(he->GetMaximumBin())<<" "<<(x_max-x_min)/2<<endl;
  TLegend *leg;
  leg= new TLegend(xleg,0.70,xleg+0.1,0.85);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);
  leg->AddEntry(he,"e-","l");
  leg->AddEntry(hmu,"#mu-","l");
  leg->AddEntry(hpi,"#pi-","l");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->Draw();
  addCaliceLogo(true);

  if(save == true){
    c->SaveAs("c_"+varname+"_"+energy+"GeV"+".eps");
    c->SaveAs("c_"+varname+"_"+energy+"GeV"+".png");
  }
}

void single_histo(TString varname, TString energy, bool save = false){

  TString filename_e = "results_folder/resolution_e-_result.root";
  TString filename_mu = "results_folder/resolution_mu-_result.root";
  TString filename_pi = "results_folder/resolution_pi-_result.root";
  
  TFile *file_e = new TFile(filename_e, "read");
  TFile *file_mu = new TFile(filename_mu, "read");
  TFile *file_pi = new TFile(filename_pi, "read");
  cout<<varname+"_e-_"+energy+"GeV"<<endl;

  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleStyle(0);
  //gStyle->SetTitleX(0.2);
  gStyle->SetMarkerSize(0.2);

  TH1 *histo_e;
  TH1 *histo_mu;
  TH1 *histo_pi;
  file_e->GetObject(varname+"_e-_"+energy+"GeV", histo_e);
  file_mu->GetObject(varname+"_mu-_"+energy+"GeV", histo_mu);
  file_pi->GetObject(varname+"_pi-_"+energy+"GeV", histo_pi);
  drawHistosTH1(varname, energy, histo_e, histo_mu, histo_pi, save);
  
  
}
