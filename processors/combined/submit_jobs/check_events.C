
void check_events(TString filename, TString detectortype, int target){

  TFile *file = new TFile(filename, "read");

  TString treename;
  if(detectortype == "LCIO2Build") treename = "ecal";
  if(detectortype == "AHCALLCIO2Build") treename = "hcal";

  TTree *tree = (TTree*)file->Get(treename);
  int entries = tree->GetEntries();
  //cout<<"File: "<<filename<<". Entries: "<<entries<<endl;  
  if( entries != target ) cout<<"ALERT! FILE "<<filename<<" has "<<entries<<" events!"<<endl;
}
