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

void findTH1values(TH1 * hA, TH1 * hB, TH1 * hC, float &x_min, float &x_max, float &y_max, bool &rescalable) {
  float max_x_A = hA->GetBinContent(hA->GetMaximumBin());
  float max_x_B = hB->GetBinContent(hB->GetMaximumBin());
  float max_x_C = hC->GetBinContent(hC->GetMaximumBin());
  
  cout<<"Max (e, mu, pi, neutron): ("<<max_x_A<<", "<<max_x_B<<", "<<max_x_C<<")"<<endl;

  int b_min_A = hA->FindFirstBinAbove(0,1);
  int b_min_B = hB->FindFirstBinAbove(0,1);
  int b_min_C = hC->FindFirstBinAbove(0,1);

  int b_max_A = hA->FindLastBinAbove(0,1);
  int b_max_B = hB->FindLastBinAbove(0,1);
  int b_max_C = hC->FindLastBinAbove(0,1);

  int b_min = min(b_min_A, min(b_min_B,b_min_C));
  int b_max = max(b_max_A, max(b_max_B,b_max_C));

  if(b_max > N_ECAL_LAYERS + 1) rescalable = true;

  x_min = hA->GetBinCenter(b_min-1);
  x_max = hA->GetBinCenter(b_max+1);

  y_max = 1.01*max(max_x_A, max(max_x_B,max_x_C));
  cout<<"(min x, max x, max y) = ("<<x_min<<", "<<x_max<<", "<<y_max<<")"<<endl;

  return 0;
}

void drawHistosTH1(TString varname, TH1 * hA, TH1 * hB, TH1 * hC, bool save = false) {
  float x_min = 0.;
  float x_max = 0.;
  float y_max = 0.;
  bool rescalable = false;
  findTH1values(hA, hB, hC, x_min, x_max, y_max, rescalable);
  hA->GetYaxis()->SetRangeUser(0,y_max);
  hA->GetYaxis()->SetMaxDigits(3);
  if(rescalable == true) hA->GetXaxis()->SetRangeUser(x_min,x_max);

  hA->SetXTitle("Value");
  hA->SetYTitle("Entries");

  auto c = new TCanvas("c_XproductionX_"+varname, "c_XproductionX_"+varname, 800, 800);
  c->cd();
  if(y_max > 1000){
    hA->SetMinimum(1);
    hB->SetMinimum(1);
    hC->SetMinimum(1);
    c->SetLogy();
  }
  hA->SetLineColor(kGray+2);
  hB->SetLineColor(kCyan-3);
  hC->SetLineColor(kRed-4);
  hA->SetLineWidth(3);
  hB->SetLineWidth(3);
  hC->SetLineWidth(3);
  
  hA->SetTitle(varname);
  hA->Draw("histo");
  hB->Draw("histosame");
  hC->Draw("histosame");

  float xleg=0.;
  if( hA->GetBinCenter(hA->GetMaximumBin()) < (x_max-x_min)/2) xleg = 0.5;
  if( hA->GetBinCenter(hA->GetMaximumBin()) > (x_max-x_min)/2) xleg = 0.2;
  TLegend *leg;
  leg= new TLegend(xleg,0.70,xleg+0.3,0.85);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->AddEntry(hA,"XsetAX","l");
  leg->AddEntry(hB,"XsetBX","l");
  leg->AddEntry(hC,"XsetCX","l");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->Draw();
  addCaliceLogo(true);

  if(save == true){
    c->SaveAs("c_XproductionX_"+varname+".eps");
    c->SaveAs("c_XproductionX_"+varname+".png");
  }
}

void single_histo_XproductionX(TString varname, bool save = false){

  TString filename_A = "../XproductionX/resolution_XresAX_result.root";
  TString filename_B = "../XproductionX/resolution_XresBX_result.root";
  TString filename_C = "../XproductionX/resolution_XresCX_result.root";

  TFile *file_A = new TFile(filename_A, "read");
  TFile *file_B = new TFile(filename_B, "read");
  TFile *file_C = new TFile(filename_C, "read");
  cout<<varname<<endl;

  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleStyle(0);
  //gStyle->SetTitleX(0.2);
  gStyle->SetMarkerSize(0.2);

  TH1 *histo_A;
  TH1 *histo_B;
  TH1 *histo_C;

  file_A->GetObject(varname+"_XpartAX", histo_A);
  file_B->GetObject(varname+"_XpartBX", histo_B);
  file_C->GetObject(varname+"_XpartCX", histo_C);

  drawHistosTH1(varname, histo_A, histo_B, histo_C, save);  
  
}
