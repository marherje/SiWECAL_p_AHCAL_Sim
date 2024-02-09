// example about structures
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
using namespace std;

map< string, vector<vector<int>> > ReadMap(TString filename) 
{
  map< string, vector<vector<int>> > map_point;
  string key[2] = {"X", "Y"};
  for(auto k : key){
    map_point[k] = vector<vector<int>>(16, vector<int>(64, -1000.));
  }

  std::ifstream reading_file(filename);
  if(!reading_file){
    cout<<" dameyo - damedame"<<endl;
  }
  
  for(int i=0; i<16; i++) {
    for(int j=0; j<64; j++) {
      map_point["X"].at(i).at(j) = -1000.;
      map_point["Y"].at(i).at(j) = -1000.;
    }
  }

  Int_t tmp_chip = 0,tmp_channel = 0;
  Float_t tmp_x0 = 0 ,tmp_y0 = 0 , tmp_x = 0 , tmp_y = 0 ;
  TString tmpst;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  while(reading_file){
    reading_file >> tmp_chip >> tmp_x0 >> tmp_y0 >> tmp_channel >> tmp_x >> tmp_y ;
    map_point["X"].at(tmp_chip).at(tmp_channel) = -tmp_x ;
    map_point["Y"].at(tmp_chip).at(tmp_channel) = -tmp_y ;
  }

  return map_point;

}


Float_t mpv_error[15][16][64];
Float_t mpv[15][16][64];


void ReadMIPCalibration(TString filename) 
{

  std::ifstream reading_file(filename);
  if(!reading_file){
    cout<<" dameyo - damedame"<<endl;
  }

  Int_t tmp_layer = 0, tmp_chip = 0,tmp_channel = 0;
  Float_t tmp_mpv = 0 ,tmp_mpv_error = 0 , tmp_ = 0 ;
  TString tmpst;
  reading_file >> tmpst >> tmpst >> tmpst ;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst>> tmpst >> tmpst >> tmpst ;
  
  while(reading_file){
    reading_file >> tmp_layer >> tmp_chip >> tmp_channel >> tmp_mpv >> tmp_mpv_error >> tmp_ >> tmp_ >> tmp_ ;
    mpv[tmp_layer][tmp_chip][tmp_channel]=tmp_mpv;
    mpv_error[tmp_layer][tmp_chip][tmp_channel]=tmp_mpv_error;
    
  }

 double ncalibrated=0.;
  double ntotal=0;
  
  for(int k=0; k<15; k++) {
    double ncalibrated_tmp=0.;
    double ntotal_tmp=0;
    for(int i=0; i<16; i++) {
      for(int j=0; j<64; j++) {
	if(mpv[k][i][j]==0 || mpv_error[k][i][j]<0) {
	  ncalibrated++;
	  ncalibrated_tmp++;
	}
	ntotal++;
	ntotal_tmp++;
      }
    }
    ncalibrated_tmp=100.-100.*ncalibrated_tmp/ntotal_tmp;
    cout<< "In Layer: "<<k <<" we read that "<<ncalibrated_tmp<<"% of channels are calibrated"<<endl;
  }
  ncalibrated=100.-100.*ncalibrated/ntotal;
  cout<< "TOTAL: In file " <<filename << " we read that "<<ncalibrated<<"% of channels are calibrated"<<endl;
  

}

Int_t masked[15][16][64];

void ReadMasked(TString filename) 
{

  std::ifstream reading_file(filename);
  if(!reading_file){
    cout<<" dameyo - damedame"<<endl;
  }

  for(int k=0; k<15; k++) {
    for(int i=0; i<16; i++) {
      for(int j=0; j<64; j++) {
	masked[k][i][j] = 0;
      }
    }
  }
  
  Int_t tmp_layer=0, tmp_chip = 0,tmp_channel = 0;
  Int_t tmp_masked = 0;
  TString tmpst;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;

  cout<<"Read Masked: "<<filename<<endl;
  
  while(reading_file){
    double masked_temp[64];
    reading_file >> tmp_layer >> tmp_chip ;
      for(int j=0; j<64; j++) reading_file >>masked_temp[j];
    for(int j=0; j<64; j++)masked[tmp_layer][tmp_chip][j]=masked_temp[j];      
  }
  
  /* string line;
  
  if(reading_file.is_open())  {
    while(!reading_file.eof())   {
      getline(reading_file, line);
      std::vector<int> masked_temp;

      size_t pos = 0;
      while ((pos = line.find(" ")) != string::npos) {
        masked_temp.push_back(std::stoi(line.substr(0, pos)));
        line.erase(0, pos + 1);
      }
      masked_temp.push_back(stoi(line.substr(0, pos)));
      masked_temp.size();
      for(int j=0; j<64; j++) masked[masked_temp.at(0)][masked_temp.at(1)][j]=masked_temp.at(j+2);  
    }
    }*/

  double nmasked=0.;
  double ntotal=0;
  for(int k=0; k<15; k++) {
    for(int i=0; i<16; i++) {
      for(int j=0; j<64; j++) {
	if(masked[k][i][j]==1) nmasked++;
	if(k==0 && i==0) cout << masked[k][i][j];
	ntotal++;
      }
    }
  }
  nmasked=100.*nmasked/ntotal;
  cout<< "In file " <<filename << " we read that "<<nmasked<<"% of channels are masked"<<endl;
  
}


//# ASU: 0 
//## ChipIndex: 0 ChipId: 0  FeedbackCap: 3 ThresholdDAC: 230 HoldDelay: 130 FSPeakTime: 2 GainSelectionThreshold: 255 
//### Ch: 0 TrigMask: 0 ChThreshold: 0 PAMask: 0

struct skiroc_t {
  int idx;
  int id;
  int n_channels=64;
  int feedback_cap;
  int threshold_dac;
  int hold_delay;
  int fs_peak_time;
  int gain_selection_threshold;
  int mask[64];
  int chn_threshold_adj[64];
  int preamplifier_mask[64];
  int default_ch_threshold=0;
  int default_trig_mask=0;
  int default_pa_mask=0;
} skiroc;

struct asu_t {
  int idx;
  int n_chips=16;
  skiroc_t skiroc[16];
} asu;

struct slab_t {
  int idx;
  int add;
  int ser_num;
  string slb_fpga;
  int nb_asus;
  asu_t asu[5];
} slab;


struct detector_t {
  string sl_soft_v;
  string date_run_unix;
  string date_run_date;
  string date_run_time;
  string usb_sernum;
  string fpga_ver;
  int external_clock;
  int n_core_daughters;
  int trigger_type;
  int acq_window_source;
  int acq_window;
  int acq_delay;
  int start_acq_delay;
  int ext_sig_level;
  int ext_sig_edge;
  string core_daughter_fpga_ver[2];
  int core_daughter_n_slabs[2];
  slab_t slab[2][15];
  bool isCore;
} detector;

bool search_string(const std::string &line, const std::string &str, std::string &par){

  std::size_t found  = line.find(str);

  if(found!=std::string::npos) {

    std::size_t found1 = line.find_first_of(":",found);
    std::size_t found2 = line.find_first_not_of(' ',found1+1);
    std::size_t found3 = line.find_first_of(' ',found2+1);

    std::string result = line.substr(found2,found3-found2);

    par = result;

    return true;

  }else{
    return false;
  }

}

bool search_string(const std::string &line, const std::string &str, int &par){

  std::size_t found  = line.find(str);

  if(found!=std::string::npos) {

    std::size_t found1 = line.find_first_of(":",found);
    std::size_t found2 = line.find_first_not_of(' ',found1+1);
    std::size_t found3 = line.find_first_of(' ',found2+1);

    std::string result = line.substr(found2,found3-found2);

    par = std::stoi(result);

    return true;

  }else{
    return false;
  }

}

bool search_string(const std::string &line, const std::string &str, bool &par){

  std::size_t found  = line.find(str);

  if(found!=std::string::npos) {

    std::size_t found1 = line.find_first_of(":",found);
    std::size_t found2 = line.find_first_not_of(' ',found1+1);
    std::size_t found3 = line.find_first_of(' ',found2+1);

    std::string result = line.substr(found2,found3-found2);

    if(result=="SL_COREMODULE_INTERFACE"){
      par = true;
    }else{
      par=false;
    }

    return true;

  }else{
    return false;
  }

}

bool search_unixT_string(const std::string &line, const std::string &str, std::string &par){

  std::size_t found  = line.find(str);

  if(found!=std::string::npos) {

    std::size_t found1 = line.find_first_of("123456789",found);
    std::size_t found2 = line.find_first_of(' ',found1);

    std::string result = line.substr(found1,found2-found1);

    par = result;

    return true;

  }else{
    return false;
  }

}

bool search_string_nocolon(const std::string &line, const std::string &str, int &par){

  std::size_t found  = line.find(str);

  if(found!=std::string::npos) {

    std::size_t found1 = line.find_first_of("0123456789",found);
    std::size_t found2 = line.find_first_of(' ',found1);

    std::string result = line.substr(found1,found2-found1);
    par = std::stoi(result);

    return true;

  }else{
    return false;
  }

}

void read_configuration_file(TString filename="Run_Settings.txt", bool debug=true) {

  std::ifstream reading_file(filename);
  if(!reading_file){
    std::cout<<" dameyo - damedame"<<std::endl;
  }

  std::string line;
  detector.isCore=false;
  int nconnasu(0);
  int core_ctr=0;

  int daughter_ctr=0;
  int slab_ctr=0;
  int asus_ctr=0;
  int chip_ctr=0;
  int ch_ctr=0;


  // TEST USE
  int tmp_count=0;


  while ( getline(reading_file, line) ) {

    if (line.c_str()[0] != 't') {

      search_string(line,"ILC SL-SOFTWARE VERSION:", detector.sl_soft_v);
      search_unixT_string(line,"UnixTime =",detector.date_run_unix);
      search_unixT_string(line,"date =",detector.date_run_date);
      search_unixT_string(line,"time =",detector.date_run_time);

      search_string(line,"SYSTEM_TYPE:", detector.isCore);

      search_string(line,"USB_SerNum:", detector.usb_sernum);
      search_string(line,"FPGA_Version:", detector.fpga_ver);
    
      if(detector.isCore){
	search_string(line,"NB_Of_Core_Daughters:", detector.n_core_daughters);
      }else{
	search_string(line,"Nb_Of_Connected_ASUs:", detector.slab[0][0].nb_asus);
      }
    
      search_string(line,"EXT_CLOCK:", detector.external_clock);

      search_string(line,"TriggerType:", detector.trigger_type);
      search_string(line,"ACQWindowSource:", detector.acq_window_source);
      search_string(line,"ACQWindow:", detector.acq_window);
      search_string(line,"DelayBetweenCycle:", detector.acq_delay);
      search_string(line,"DelayForStartAcq:", detector.start_acq_delay);
      search_string(line,"ExtSigLevel:", detector.ext_sig_level);
      bool isLast = search_string(line,"ExtSigEdge:", detector.ext_sig_edge);

      if(isLast) break;

    } // initial check

  } // loop line


  if(debug){

    cout << "SL-SOFTWARE VERSION: " << detector.sl_soft_v << "\n";
    cout << "UnixTime = " << detector.date_run_unix <<", date = "<< detector.date_run_date <<", time = "<< detector.date_run_time <<"\n";

    cout << "SYSTEM_TYPE: " << detector.isCore << "\n";
    cout << "USB_SerNum: " << detector.usb_sernum << "\n";
    cout << "FPGA_Version: " << detector.fpga_ver << "\n";
    cout << "NB_Of_Core_Daughters: " << detector.n_core_daughters << "\n";
    cout << "EXT_CLOCK: " << detector.external_clock << "\n";

    cout << "TriggerType: " << detector.trigger_type << "\n";
    cout << "ACQWindowSource: " << detector.acq_window_source << "\n";
    cout << "ACQWindow: " << detector.acq_window << "\n";
    cout << "DelayBetweenCycle: " << detector.acq_delay << "\n";
    cout << "DelayForStartAcq: " << detector.start_acq_delay << "\n";
    cout << "ExtSigLevel: " << detector.ext_sig_level << "\n";
    cout << "ExtSigEdge: " << detector.ext_sig_edge << "\n";

  }

  
  if(detector.isCore) {

    for (int i = 0; i < detector.n_core_daughters; ++i)
    {
      int icore=-1;
      getline(reading_file, line);
      search_string(line,"CORE_DAUGHTER:",icore);
      search_string(line,"FPGA_Version:",detector.core_daughter_fpga_ver[i]);
      search_string(line,"Nb_Of_Connected_Slabs:",detector.core_daughter_n_slabs[i]);
    }

    if(debug){

      for(int i=0; i < detector.n_core_daughters; ++i){
        cout << "FPGA_Version: " << detector.core_daughter_fpga_ver[i] << endl;
        cout << "Nb_Of_Connected_Slabs: " << detector.core_daughter_n_slabs[i] << endl;
      }

    }

  }else{
    
    detector.n_core_daughters         = 1;
    detector.core_daughter_n_slabs[0] = 1;

  }


  for (int i = 0; i < detector.n_core_daughters; i++) {
    
    for (int j = 0; j < detector.core_daughter_n_slabs[i]; j++) {

      if (detector.isCore) {
        
        getline(reading_file, line);
        search_string(line,"SlabIdx:",detector.slab[i][j].idx);
        search_string(line,"SlabAdd:",detector.slab[i][j].add);
        search_string(line,"SL_Board_SerNum:",detector.slab[i][j].ser_num);
        search_string(line,"FPGA_Version:",detector.slab[i][j].slb_fpga);
        search_string(line,"Nb_Of_Connected_ASUs:",detector.slab[i][j].nb_asus);

        if (debug) cout << " Daughter: " << i << " SlabIdx: " << detector.slab[i][j].idx << " SlabAdd: " << detector.slab[i][j].add << " SL_Board_SerNum: " << detector.slab[i][j].ser_num << " FPGA_Version: " << detector.slab[i][j].slb_fpga << " Nb_Of_Connected_ASUs: " << detector.slab[i][j].nb_asus << endl;
      }

      // Calib Pulse Parameter
      getline(reading_file,line);

      //CRP
      for (int k = 0; k < detector.slab[i][j].nb_asus; k++) {

        getline(reading_file,line);
        search_string(line,"ASU:",detector.slab[i][j].asu[k].idx);

        if (debug) cout << "# ASU: " << detector.slab[i][j].asu[k].idx << " " << endl;

        for (int ichip = 0; ichip < detector.slab[i][j].asu[k].n_chips; ichip++) {

          getline(reading_file,line);

          search_string(line,"ChipIndex:",detector.slab[i][j].asu[k].skiroc[ichip].idx);
          search_string(line,"ChipId:",detector.slab[i][j].asu[k].skiroc[ichip].id);
          search_string(line,"FeedbackCap:",detector.slab[i][j].asu[k].skiroc[ichip].feedback_cap);
          search_string(line,"ThresholdDAC:",detector.slab[i][j].asu[k].skiroc[ichip].threshold_dac);
          search_string(line,"HoldDelay:",detector.slab[i][j].asu[k].skiroc[ichip].hold_delay);
          search_string(line,"FSPeakTime:",detector.slab[i][j].asu[k].skiroc[ichip].fs_peak_time);
          search_string(line,"GainSelectionThreshold:",detector.slab[i][j].asu[k].skiroc[ichip].gain_selection_threshold);
          
          if (debug) cout << "## ChipIndex: " << detector.slab[i][j].asu[k].skiroc[ichip].idx << " ChipId: " << detector.slab[i][j].asu[k].skiroc[ichip].id << "  FeedbackCap: " << detector.slab[i][j].asu[k].skiroc[ichip].feedback_cap << " ThresholdDAC: " << detector.slab[i][j].asu[k].skiroc[ichip].threshold_dac << " HoldDelay: " << detector.slab[i][j].asu[k].skiroc[ichip].hold_delay << " FSPeakTime: " << detector.slab[i][j].asu[k].skiroc[ichip].fs_peak_time << " GainSelectionThreshold: " << detector.slab[i][j].asu[k].skiroc[ichip].gain_selection_threshold << endl;

          for (int ichn = 0; ichn < detector.slab[i][j].asu[k].skiroc[ichip].n_channels; ichn++) {

            getline(reading_file,line);
            search_string(line,"TrigMask:",detector.slab[i][j].asu[k].skiroc[ichip].mask[ichn]);
            search_string(line,"ChThreshold:",detector.slab[i][j].asu[k].skiroc[ichip].chn_threshold_adj[ichn]);
            search_string(line,"PAMask:",detector.slab[i][j].asu[k].skiroc[ichip].preamplifier_mask[ichn]);

            if (debug) cout << "### Ch: " << ichn << " TrigMask: " << detector.slab[i][j].asu[k].skiroc[ichip].mask[ichn] << " ChThreshold: " << detector.slab[i][j].asu[k].skiroc[ichip].chn_threshold_adj[ichn] << " PAMask: " << detector.slab[i][j].asu[k].skiroc[ichip].preamplifier_mask[ichn] << endl;

          } //End of loop over channels
        
        } //End of loop over chips of an ASU
      
      } //End of loop over ASUs
    
    } //End of loop over slabs connected to a core daughter

  } //End of loop over core daughters

}

void mask_full_chip(int idaughter, int islab, int iasu, int ichip) {
  for(int i=0; i<detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].n_channels; i++) detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].mask[i]=1;
}

void disable_trig_otherrows(int idaughter, int islab, int iasu, int ichip, int jrow) {
  //row 0 --> chns 0-7
  //row 1--> chns 8-15 
  //etc 

  for(int irow=0; irow<8; irow++) {
    if(jrow==irow) continue;
    cout<<"Disable trigger for rows "<< irow<<endl;
    for(int i=irow*8; i<(irow*8+8); i++) {
      cout<<idaughter<<" "<<islab<<" "<<iasu<<" "<< ichip<<" "<<i <<endl;
      detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].mask[i]=1;
    }
  }
}

void enable_trig_row(int idaughter, int islab, int iasu, int ichip, int irow) {
  //row 0 --> chns 0-7
  //row 2--> chns 8-15
  //etc
  cout<<irow<<endl;
  for(int i=irow*8; i<(irow*8+8); i++) {
    cout<<idaughter<<" "<<islab<<" "<<iasu<<" "<< ichip<<" "<<i <<endl;
    detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].mask[i]=0;
  }
}

TString trigger_type_string(int i) {

  TString result="SELF_TRIGGER";
  if(i==0) result="FORCE_TRIGGER";
  if(i==1) result ="SELF_TRIGGER";

  return result;
}

TString external_clock_string(int i) {

  TString result="NO";
  if(i==1) result="YES";
  return result;
}

TString acq_window_source_string(int i) {
  
  //('0' = AUTO, '1'= SOFT_CMD, '2' = EXT_SIG)
  TString result="AUTO";
  if(i==0) result="AUTO";
  if(i==1) result ="SOFT_CMD";
  if(i==2) result ="EXT_SIG";

  return result;
}

TString ext_sig_level_string(int i) {
  
  //  ('0' = TTL, '1' = NIM, '-1' = NA)
  TString result="TTL";
  //  if(i==0) result="TTL";
  //if(i==1) result ="NIM";
  //if(i==-1) result ="NA";

  return result;
}


TString ext_sig_edge_string(int i) {
  
  //  ('0' = RISING_EDGE, '1' = FALLING_EDGE, '-1' = NA)
  TString result="POS";
  // if(i==0) result="RISING_EDGE";
  //if(i==1) result ="FALLING_EDGE";
  //if(i==-1) result ="NA";
  
  return result;
}


void write_configuration_file(TString filename="Run_Settings_2.txt") {

  //===== Daughter: 0 SlabIdx: 0 SlabAdd: 1 SL_Board_SerNum: -1 FPGA_Version: V2.3.7  Nb_Of_Connected_ASUs: 1 =====
  ofstream fout(filename,ios::out);

//fout<<"== SETTINGS FILE SAVED WITH ILC SL-SOFTWARE VERSION: "<<detector.sl_soft_v<<"  == DATE OF RUN: UnixTime = "<<detector.date_run_unix<<" date = "<<detector.date_run_date<<" time = "<<detector.date_run_time<<"  ==="<<endl;
//fout<<"== SYSTEM_TYPE: SL_COREMODULE_INTERFACE USB_SerNum: "<<detector.usb_sernum<<" FPGA_Version: "<<detector.fpga_ver<<"  NB_Of_Core_Daughters: "<<detector.n_core_daughters<<" =="<<endl;
//fout<<"== TriggerType: "<<detector.trigger_type<<"  ('0' = FORCE_TRIGGER, '1' = SELF_TRIGGER) ACQWindowSource: "<<detector.acq_window_source <<" ('0' = AUTO, '1'= SOFT_CMD, '2' = EXT_SIG) ACQWindow: "<< detector.acq_window <<" (ms) DelayBetweenCycle: "<< detector.acq_delay <<" (ms) ExtSigLevel: "<< detector.ext_sig_level  <<" ('0' = TTL, '1' = NIM, '-1' = NA) ExtSigEdge: "<< detector.ext_sig_edge <<" ('0' = RISING_EDGE, '1' = FALLING_EDGE, '-1' = NA) =="<<endl;

  fout<<"# AcqParams TrigType ' "<<trigger_type_string(detector.trigger_type)<<" ' AcqWindowSource ' "<<acq_window_source_string(detector.acq_window_source) <<" ' ACQWindow "<< detector.acq_window <<" DelayBetweenCycles "<< detector.acq_delay <<" DelayForStartAcq 0 ExtSigLvel ' "<< ext_sig_level_string(detector.ext_sig_level)  <<" ' ExtSigEdge ' "<< ext_sig_edge_string(detector.ext_sig_edge) <<" '"<<endl;

 //CRP 14/12/20 This is a terrible hack!!!  
 //It is to make use of the current organisation of the code
 //This clearly has to be rewritten soon the more since the code that fills the detector struct also tacitly supposed only
 //one CORE Daughter and will thus collapse as sonn as a second CORE Daughter will be added!!!!
  if(!detector.isCore) {
   std::cout << "No CORE Module present" << std::endl;
   detector.n_core_daughters = 1;
 }
 //
 //----------------------
  //=== CORE_DAUGHTER: 0 FPGA_Version: V1.3.1 Nb_Of_Connected_Slabs: 6 ===
 for(int i=0; i<detector.n_core_daughters; i++) {
     //  fout<<"=== CORE_DAUGHTER: "<<i<<" FPGA_Version: "<<detector.core_daughter_fpga_ver[i]<<" Nb_Of_Connected_Slabs: "<<detector.core_daughter_n_slabs[i]<<" ==="<<endl;
   std::cout << "In writing loop" << std::endl;
   for(int j=0; j<detector.core_daughter_n_slabs[i]; j++) {
     std::cout << "In writing loop 2" << std::endl;
     for(int k=0; k<detector.slab[i][j].nb_asus; k++) {
       std::cout << "In writing loop 3" << std::endl;
	 //	fout<<"# SlabSerNum "<<detector.slab[i][j].ser_num<<" SlabAdd "<< detector.slab[i][j].add <<" Asu "<<detector.slab[i][j].asu[k].idx<<endl;
       fout<<"# SlabSerNum 0.0 SlabAdd "<< detector.slab[i][j].add <<" Asu "<<detector.slab[i][j].asu[k].idx<<endl;
       for(int ichip=0; ichip<detector.slab[i][j].asu[k].n_chips; ichip++) {
        fout<<"## ChipIndex "<< detector.slab[i][j].asu[k].skiroc[ichip].idx<<" FeedbackCap "<<detector.slab[i][j].asu[k].skiroc[ichip].feedback_cap<<" ThresholdDAC "<<detector.slab[i][j].asu[k].skiroc[ichip].threshold_dac <<" HoldDelay "<<  detector.slab[i][j].asu[k].skiroc[ichip].hold_delay <<" FSPeakTime "<< detector.slab[i][j].asu[k].skiroc[ichip].fs_peak_time <<" GainSelectionThreshold "<< detector.slab[i][j].asu[k].skiroc[ichip].gain_selection_threshold <<" DefaultChThreshold "<<detector.slab[i][j].asu[k].skiroc[ichip].default_ch_threshold <<" DefaultTrigMask "<<detector.slab[i][j].asu[k].skiroc[ichip].default_trig_mask<<" DefaultPAMask "<< detector.slab[i][j].asu[k].skiroc[ichip].default_pa_mask <<endl;

        for(int ichn=0; ichn<detector.slab[i][j].asu[k].skiroc[ichip].n_channels; ichn++) {
          fout<<"### Ch "<<ichn<<" TrigMask "<<detector.slab[i][j].asu[k].skiroc[ichip].mask[ichn]<<" ChThreshold "<< detector.slab[i][j].asu[k].skiroc[ichip].chn_threshold_adj[ichn] <<" PAMask "<< detector.slab[i][j].asu[k].skiroc[ichip].preamplifier_mask[ichn] <<" "<<endl;
        }
      }
    }
     }//loop over slabs
   }//loop over core daughters

 }

void disable_PA_mask(){

  for(int idaughter=0; idaughter < detector.n_core_daughters; idaughter++) {

    for(int islab=0; islab<detector.core_daughter_n_slabs[idaughter]; islab++) {

      for (int iasu = 0; iasu < detector.slab[idaughter][islab].nb_asus; ++iasu) {

        for (int ichip = 0; ichip < detector.slab[idaughter][islab].asu[iasu].n_chips; ++ichip) {

          for(int ichn=0; ichn<detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].n_channels; ichn++) {

            if( detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].mask[ichn]==1 ) {
              detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].preamplifier_mask[ichn] = 1;
            }

          } // channel

        }   // chip

      }     // asu

    }       // slab

  }         // daughter

}
