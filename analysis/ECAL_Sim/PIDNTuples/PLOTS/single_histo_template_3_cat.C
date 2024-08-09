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
  hA->GetYaxis()->SetTitleOffset(1.3);
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
  hA->DrawNormalized("histo");
  hB->DrawNormalized("histosame");
  hC->DrawNormalized("histosame");

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

void single_histo_3_temp(TString varname, bool save = false){

  TString filename_A = "../XproductionX/resolution_XresAX_result.root";
  TString filename_B = "../XproductionX/resolution_XresBX_result.root";
  TString filename_C = "../XproductionX/resolution_XresCX_result.root";

  TFile *file_A = new TFile(filename_A, "read");
  TFile *file_B = new TFile(filename_B, "read");
  TFile *file_C = new TFile(filename_C, "read");
  cout<<varname<<endl;

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetMarkerSize(0.2);

  TTree *tree_A = (TTree*)file_A->Get("ntp");
  TTree *tree_B = (TTree*)file_B->Get("ntp");
  TTree *tree_C = (TTree*)file_C->Get("ntp");

  TBranch *branch_A = 0;
  TBranch *branch_B = 0;
  TBranch *branch_C = 0;

  float value_A = 0.;
  float value_B = 0.;
  float value_C = 0.;

  tree_A->SetBranchAddress(varname, &value_A, &branch_A);
  tree_B->SetBranchAddress(varname, &value_B, &branch_B);
  tree_C->SetBranchAddress(varname, &value_C, &branch_C);

  int entries_A = tree_A->GetEntries();
  int entries_B = tree_B->GetEntries();
  int entries_C = tree_C->GetEntries();
  cout<<"entries: "<<entries_A<<" "<<entries_B<<" "<<entries_C<<endl;

  vector<float> allvalues;
  for (int i = 0; i < entries_A; i++) {
    tree_A->GetEntry(i);
    allvalues.push_back(value_A);
  }
  for (int i = 0; i < entries_B; i++) {
    tree_B->GetEntry(i);
    allvalues.push_back(value_B);
  }
  for (int i = 0; i < entries_C; i++) {
    tree_C->GetEntry(i);
    allvalues.push_back(value_C);
  }
  float xmin = *std::min_element(allvalues.begin(), allvalues.end());;
  float xmax = *std::max_element(allvalues.begin(), allvalues.end());;
  cout<<"min: "<<xmin<<", max: "<<xmax<<endl;

  TH1F *histo_A = new TH1F("hA", "hA", 100, xmin, xmax);
  TH1F *histo_B = new TH1F("hB", "hB", 100, xmin, xmax);
  TH1F *histo_C = new TH1F("hC", "hC", 100, xmin, xmax);

  for (int i = 0; i < entries_A; i++) {
    tree_A->GetEntry(i);
    //cout<<value_A<<endl;
    histo_A->Fill(value_A);
  }

  for (Int_t i=0; i<entries_B; i++){
    tree_B->GetEntry(i);
    //cout<<value_B<<endl;
    histo_B->Fill(value_B);
  }

  for (Int_t i=0; i<entries_C; i++){
    tree_C->GetEntry(i);
    //cout<<value_C<<endl;
    histo_C->Fill(value_C);
  }

  drawHistosTH1(varname, histo_A, histo_B, histo_C, save);

}
