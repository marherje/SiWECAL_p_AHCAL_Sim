#define N_ECAL_LAYERS 15

void addCaliceLogo(bool WIP = true){
  TImage *img = TImage::Open("../../style/CALICELogo_18pc.png");
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
  int b_max = max(b_max_e, min(b_max_mu,b_max_pi));

  if(b_max > N_ECAL_LAYERS + 1) rescalable = true;

  x_min = he->GetBinCenter(b_min-1);
  x_max = he->GetBinCenter(b_max+1);

  y_max = 1.01*max(max_x_e, max(max_x_mu,max_x_pi));
  cout<<"(min x, max x, max y) = ("<<x_min<<", "<<x_max<<", "<<y_max<<")"<<endl;

  return 0;
}

void drawHistosTH1(TString varname, TH1 * he, TH1 * hmu, TH1 * hpi, bool save = false) {
  float x_min = 0.;
  float x_max = 0.;
  float y_max = 0.;
  bool rescalable = false;
  findTH1values(he, hmu, hpi, x_min, x_max, y_max, rescalable);
  he->GetYaxis()->SetRangeUser(0,y_max);
  he->GetYaxis()->SetMaxDigits(3);
  if(rescalable == true) he->GetXaxis()->SetRangeUser(x_min,x_max);

  he->SetXTitle("Value");
  he->SetYTitle("Entries");

  auto c = new TCanvas("c_PID_6_60_150_GeV_e_"+varname, "c_PID_6_60_150_GeV_e_"+varname, 800, 800);
  c->cd();
  he->SetLineColor(kGray+2);
  hmu->SetLineColor(kCyan-3);
  hpi->SetLineColor(kRed-4);
  he->SetLineWidth(3);
  hmu->SetLineWidth(3);
  hpi->SetLineWidth(3);

  he->SetTitle(varname);
  he->Draw("histo");
  hmu->Draw("histosame");
  hpi->Draw("histosame");

  float xleg=0.575;
  if( he->GetMaximumBin() < he->GetNbinsX()/2) xleg = 0.65;
  TLegend *leg;
  leg= new TLegend(xleg,0.70,xleg+0.2,0.85);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);
  leg->AddEntry(he,"e- (6 GeV)","l");
  leg->AddEntry(hpi,"e- (60 GeV)","l");
  leg->AddEntry(hmu,"e- (150 GeV)","l");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->Draw();
  addCaliceLogo(true);

  if(save == true){
    c->SaveAs("c_PID_6_60_150_GeV_e_"+varname+".eps");
    c->SaveAs("c_PID_6_60_150_GeV_e_"+varname+".png");
  }
}

void single_histo_PID_6_60_150_GeV_e(TString varname, bool save = false){

  TString filename_e = "../PID_6_60_150_GeV_e/resolution_6_e-_result.root";
  TString filename_pi = "../PID_6_60_150_GeV_e/resolution_60_e-_result.root";
  TString filename_mu = "../PID_6_60_150_GeV_e/resolution_150_e-_result.root";
  
  TFile *file_e = new TFile(filename_e, "read");
  TFile *file_pi = new TFile(filename_pi, "read");
  TFile *file_mu = new TFile(filename_mu, "read");
  cout<<varname<<endl;

  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleStyle(0);
  //gStyle->SetTitleX(0.2);
  gStyle->SetMarkerSize(0.2);

  TH1 *histo_e;
  TH1 *histo_pi;
  TH1 *histo_mu;

  file_e->GetObject(varname+"_e-", histo_e);
  file_pi->GetObject(varname+"_e-", histo_pi);
  file_mu->GetObject(varname+"_e-", histo_mu);

  drawHistosTH1(varname, histo_e, histo_mu, histo_pi, save);  
  
}
