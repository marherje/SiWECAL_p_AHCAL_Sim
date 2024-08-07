#include "LCIO2BuildProcessor.hh"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <map>

// ---- LCIO Headers
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include "IMPL/CalorimeterHitImpl.h"
#include "IMPL/SimCalorimeterHitImpl.h"
#include <IMPL/LCFlagImpl.h>
#include "UTIL/LCTypedVector.h"
#include "UTIL/BitField64.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#endif // MARLIN_USE_AIDA

//ROOT
#include <TH1.h>
#include <TH2.h>
#include "TF1.h"
#include <TString.h>
#include <TFile.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TSpectrum.h>
#include "TMath.h"
#include "TROOT.h"
#include <TTree.h>
#include <TInterpreter.h>


#include <fstream>
#include <sstream>

//  //From Adrian
#include <vector>
#include <algorithm>


using namespace lcio ;
using namespace marlin ;
using namespace std;

namespace CALICE
{

  LCIO2BuildProcessor aLCIO2BuildProcessor;

  LCIO2BuildProcessor::LCIO2BuildProcessor() : Processor("LCIO2BuildProcessor")
  {

    vector<string> calorimInpCollectionsExample;
    registerProcessorParameter("Input_Collections",
                               "vector of collections for input",
                               _calorimInpCollections,
                               calorimInpCollectionsExample);

    string outFileNameExample;
    registerProcessorParameter("OutputBuildFile",
                               "Name to identify output file",
                               _outFileName,
                               outFileNameExample);
    
    vector<string> siThicknessesExample;
    registerProcessorParameter("SiThicknesses",
                               "vector of Silicon thicknesses per layer",
                               _siThicknesses,
                               siThicknessesExample);
    
    vector<string> FixedPosZExample;
    registerProcessorParameter("FixedPosZ",
                               "vector of z layer positions",
                               _FixedPosZ,
                               FixedPosZExample);
    
    registerProcessorParameter("NSlabs",
                               "Number of slab",
                               _NSlabs,
                               int(9999));

    registerProcessorParameter("FirstSlabPosZ",
                               "Position of first slab",
                               _FirstSlabPosZ,
                               float(9999));

    registerProcessorParameter("SlabSpacing",
                               "Separation between slabs",
                               _SlabSpacing,
                               float(9999));
    
    vector<string> GeV2MIPFactorsExample;
    registerProcessorParameter("GeV2MIPFactors",
                               "vector of MIP conversion factors",
                               _GeV2MIPFactors,
                               GeV2MIPFactorsExample);
    
    registerProcessorParameter("ConversionGeV2MIP",
                               "Do conversion to MIP",
                               _ConversionGeV2MIP,
                               false);
    
    registerProcessorParameter("MIPThreshold",
                               "Minimum energy (in MIPs) to count a hit",
                               _MIPThreshold,
                               float(0.));

    // vector<string> MapFilenamesExample = {"/home/llr/ilc/jimenez/Projects/Simulations/SiWECAL-Sim/processors/LCIO2BuildProcessor/mapping/fev10_chip_channel_x_y_mapping.txt", "/home/llr/ilc/jimenez/Projects/Simulations/SiWECAL-Sim/processors/LCIO2BuildProcessor/mapping/fev11_cob_rotate_chip_channel_x_y_mapping.txt"};
    vector<string> MapFilenamesExample = {"", ""}; 
      //{"/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/mapping/fev10_chip_channel_x_y_mapping.txt", "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/mapping/fev11_cob_rotate_chip_channel_x_y_mapping.txt"};
    registerProcessorParameter("MappingFiles",
                               "Files mapping hit position with cell chan and chip",
                               _MapFilenames,
                               MapFilenamesExample);
    
    // vector<int> SlabMapIndicesExample = {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1};
    vector<int> SlabMapIndicesExample = {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0};
    registerProcessorParameter("SlabMapIndices",
                               "vector indices of maps (from MappingFiles parameters) to be used per slab",
                               _SlabMapIndices,
                               SlabMapIndicesExample);    


    float HalfCenterGapExample = 1.05; // This dead space should be included in the generation (3.8 - 5.5/2)
    registerProcessorParameter("HalfCenterGap",
                               "Half size of gap between wafers",
                               _HalfCenterGap,
                               HalfCenterGapExample);
    
    // registerProcessorParameter("hitType",
    //                            "Hit type (SimCalorimeterHit or CalorimeterHit)",
    //                            _hitType,
    //                            "SimCalorimeterHit");
    
  }

  /************************************************************************************/

  void LCIO2BuildProcessor::init()
  {
    streamlog_out(DEBUG) << "init called" << std::endl ;

    printParameters();

    _nRun = 0 ;
    _nEvt = 0 ;
    _rootout = new TFile(_outFileName.c_str(), "RECREATE");
    _treeout = new TTree("ecal", "From SLCIO");
    _treeout->Branch("event", &event);
    _treeout->Branch("spill", &spill);
    _treeout->Branch("cycle", &cycle);
    _treeout->Branch("bcid", &bcid);
    _treeout->Branch("bcid_first_sca_full", &bcid_first_sca_full);
    _treeout->Branch("bcid_merge_end", &bcid_merge_end);
    _treeout->Branch("id_run", &id_run);
    _treeout->Branch("id_dat", &id_dat);
    _treeout->Branch("nhit_slab", &nhit_slab);
    _treeout->Branch("nhit_chip", &nhit_chip);
    _treeout->Branch("nhit_chan", &nhit_chan);
    _treeout->Branch("nhit_len", &nhit_len);

    // _treeout->Branch("sum_hg", &sum_hg);
    _treeout->Branch("sum_energy", &sum_energy);
    _treeout->Branch("sum_energy_og", &sum_energy_og);

    _treeout->Branch("hit_slab", &hit_slab);
    _treeout->Branch("hit_chip", &hit_chip);
    _treeout->Branch("hit_chan", &hit_chan);
    _treeout->Branch("hit_sca", &hit_sca);
    _treeout->Branch("hit_adc_high", &hit_adc_high);
    _treeout->Branch("hit_adc_low", &hit_adc_low);
    _treeout->Branch("hit_n_scas_filled", &hit_n_scas_filled);
    _treeout->Branch("hit_isHit", &hit_isHit);
    _treeout->Branch("hit_isMasked", &hit_isMasked);
    _treeout->Branch("hit_isCommissioned", &hit_isCommissioned);

    _treeout->Branch("hit_energy", &hit_energy);
    _treeout->Branch("hit_energy_w", &hit_energy_w);
    _treeout->Branch("hit_energy_og", &hit_energy_og);
    _treeout->Branch("hit_x", &hit_x);
    _treeout->Branch("hit_y", &hit_y);
    _treeout->Branch("hit_z", &hit_z);

    
    if(_FixedPosZ.size() != 0) {
      for (int ilayer = 0; ilayer < int(_FixedPosZ.size()); ilayer++){
        _FixedPosZ_float.push_back(stof(_FixedPosZ[ilayer]));
      }
    }
    else {
      if(_NSlabs == 9999){ streamlog_out(ERROR)<<"Missing number of slabs"<<std::endl;}
      else {
        if(_FirstSlabPosZ == 9999.){ streamlog_out(ERROR)<<"Missing first slab position"<<std::endl;}
        else if(_SlabSpacing == 9999.){ streamlog_out(ERROR)<<"Missing slab separation"<<std::endl;}
        else{
          for (int ilayer = 0; ilayer < _NSlabs; ilayer++){
            _FixedPosZ_float.push_back(_FirstSlabPosZ + _SlabSpacing*ilayer);
          }
        }
      }
    }

    for (int ilayer = 0; ilayer < int(_FixedPosZ.size()); ilayer++){
      if( _siThicknesses.size() != 1 ) {
        if(_FixedPosZ.size() != _siThicknesses.size()) {streamlog_out(ERROR)<< "Si Thicknesses vector size don't match the number of layers "<<std::endl;}
        else {
          _siThicknesses_float.push_back(0.001*stof(_siThicknesses[ilayer])); // Converted from microns to mm
        }
      }
      else {
        try{
          if(_siThicknesses.size() == 1) _siThicknesses_float.push_back(0.001*stof(_siThicknesses[0])); //Converted from microns to mm
          else throw 404;
        }
        catch(...){streamlog_out(ERROR)<< "Missing Si Thicknesses info "<<std::endl;
        }
      }
    }

    if(_ConversionGeV2MIP==true){
      for (int ilayer = 0; ilayer < int(_FixedPosZ.size()); ilayer++){
	if( _GeV2MIPFactors.size() > 1 ) {
	  if(_FixedPosZ.size() != _GeV2MIPFactors.size()) {streamlog_out(ERROR)<< "MIP vector size don't match the number of layers "<<std::endl;}
	  else {
	    _GeV2MIPFactors_float.push_back(stof(_GeV2MIPFactors[ilayer]));
	  }
	}
	else {
	  try{
	    if(_GeV2MIPFactors.size() == 1) _GeV2MIPFactors_float.push_back(stof(_GeV2MIPFactors.at(0)));
	    else throw 404;
	  }
	  catch(...){streamlog_out(ERROR)<< "Missing MIP info "<<std::endl;
	  }
	}
      }
    }
    
    streamlog_out(DEBUG)<<"PARAMETER READING"<<std::endl;
    streamlog_out(DEBUG)<<"GeV2MIP size: "<<_GeV2MIPFactors.size()<<", GeV2MIP: "<<_GeV2MIP<<std::endl;
    streamlog_out(DEBUG)<<"Threshold energy (MIPS): "<<_MIPThreshold<<std::endl;
    streamlog_out(DEBUG)<<"FixedPosZ size: "<<_FixedPosZ.size()<<std::endl;
    streamlog_out(DEBUG)<<"siThicknesses size: "<<_siThicknesses.size()<<std::endl;
    streamlog_out(DEBUG)<<"N Slabs: "<<_NSlabs<<", First slab pos Z: "<<_FirstSlabPosZ<<", Slab spacing: "<<_SlabSpacing<<std::endl;

    
    // _FixedPosZ = {6.225,  21.225,  36.15,  51.15,  66.06,  81.06,  96.06, \
    //               111.15, 126.15, 141.15, 156.15, 171.06, 186.06, 201.06,\
    //               216.06};
    // ILD mode
    // _FixedPosZ = {6.225,  21.225,  36.15,  51.15,  66.06,  81.06,  96.06,\
    //               111.15, 126.15, 141.15, 156.15, 171.06, 186.06, 201.06,\
    //               216.06, 216., 231., 246., 261., 276., 291., 306., 321., 336.};
    _printType = true;

    // Mapping (should be done with nice dict...)
    // vector<vector<vector<float>>> _maps(_MapFilenames.size(), vector<vector<float>>(1024, vector<float>(2)));
    // Initializing to zeros, needs to be done in a smarter way
    vector<float> zeros = {0., 0.};
    vector<vector<float>> this_map;
    for (int imap = 0; imap < _MapFilenames.size(); imap++) {
      for (int icell = 0; icell < 1024; icell++) this_map.push_back(zeros);
      _maps.push_back(this_map); 
    }
    fstream map_file;
    string space_delimiter = " ";
    size_t pos;
    vector<string> words;

    for (int imap = 0; imap < _MapFilenames.size(); imap++) {
      map_file.open(_MapFilenames[imap], ios::in);
      vector<vector<float>> this_map;
      if (map_file.is_open()){   
        string line;
        getline(map_file, line); //skip the first line
        int iline = 0;
        while(getline(map_file, line)){ 
          auto start = 0U;
          auto end = line.find(space_delimiter);
          while (end != std::string::npos) {
              words.push_back(line.substr(start, end - start));
              start = end + space_delimiter.length();
              end = line.find(space_delimiter, start);
          }
          words.push_back(line.substr(start));
          // There's a sign flip in coordinates in mapping (hence the -1*)
          _maps[imap][64 * stoi(words[0]) + stoi(words[3])][0] = -1 * stof(words[4]);
          _maps[imap][64 * stoi(words[0]) + stoi(words[3])][1] = -1 * stof(words[5]);
          words.clear();
        }

        map_file.close();
      }
      else streamlog_out(WARNING)<<"Mapping file not found! Add input or check src and mapping folders"<<endl;
    }
  }


  /************************************************************************************/
  void LCIO2BuildProcessor::processRunHeader( LCRunHeader* run)
  {
    _nRun++ ;
  }

  void LCIO2BuildProcessor::processEvent( LCEvent * evt )
  {
    
    int evtNumber = evt->getEventNumber();
    if ((evtNumber % 1000) == 0)
    streamlog_out(DEBUG) << " \n ---------> Event: "<<evtNumber<<"!Collections!! <-------------\n" << std::endl;

    //Get the list of collections in the event
    const std::vector<std::string> *cnames = evt->getCollectionNames();
    
    std::vector<float> this_energyCont;

    // int this_nMC;
    for(unsigned int icol = 0; icol < _calorimInpCollections.size(); icol++)
    {
      if( std::find(cnames->begin(), cnames->end(), _calorimInpCollections.at(icol)) != cnames->end() )
      {
        string _inputColName = _calorimInpCollections.at(icol);

        LCCollection *inputCalorimCollection = 0;
        try
        {
          inputCalorimCollection = evt->getCollection(_inputColName);
        }
        catch (EVENT::DataNotAvailableException &e)
        {
          streamlog_out(WARNING)<< "missing collection "
          <<_inputColName<<endl<<e.what()<<endl;
        }
        if (inputCalorimCollection != 0){

          int noHits = inputCalorimCollection->getNumberOfElements();

          event = evtNumber;
          spill = -1;
          cycle = -1;
          bcid = -1;
          bcid_first_sca_full = -999;
          bcid_merge_end = -1;
          id_run = -1;
          id_dat = -1;
          nhit_len = noHits;

          sum_energy_og = 0.;
          sum_energy = 0.;

          vector<int> slabs_hit, chans_hit, chips_hit;
          int i_slab;
          float gap_hit_x, gap_hit_y;

          // Hits loop
          for (int i = 0; i < noHits; i++)
          {
            // // Following line works, the pushbacks
            SimCalorimeterHit *aHit = dynamic_cast<SimCalorimeterHit*>(inputCalorimCollection->getElementAt(i));
            // //auto *aHit = hitCast(inputCalorimCollection->getElementAt(i));
            // // hitCast(*aHit);

            // This is in DigiLCIO2Build
            // CalorimeterHit *aHit = dynamic_cast<CalorimeterHit*>(inputCalorimCollection->getElementAt(i));

            hit_sca.push_back(-1);
            hit_adc_high.push_back(-1);
            hit_adc_low.push_back(-1);
            hit_n_scas_filled.push_back(1);
            hit_isHit.push_back(1);
            hit_isMasked.push_back(0);
            hit_isCommissioned.push_back(1);
                
	    //streamlog_out(DEBUG)<<"CellID0: "<<aHit->getCellID0()<<", CellID1: "<<aHit->getCellID1()<<endl;

	    // Insertion Jesus                                                                  
	    vector<float> hitZtoslab;
	    float smallestdistance=9999.;
	    float hitZ=aHit->getPosition()[2];
	    
	    if(_FixedPosZ.size() == 0) streamlog_out(ERROR)<<"Missing slabs info"<<std::endl;

	    for (int islab = 0; islab < _FixedPosZ_float.size(); islab++){
	      hitZtoslab.push_back(abs(_FixedPosZ_float[islab]-hitZ));
	    }
	    for (int islab = 0; islab < _FixedPosZ_float.size(); islab++){
	      streamlog_out(DEBUG)<<"layer "<<islab<<" distance: "<<hitZtoslab.at(islab)<<endl;
	      if ( (_FixedPosZ_float[islab] > (aHit->getPosition()[2] - _siThicknesses_float.at(islab)) ) &&
		   (_FixedPosZ_float[islab] < (aHit->getPosition()[2] + _siThicknesses_float.at(islab)) ) ) {
                i_slab = islab;
		break;
	      }
              else if (_FixedPosZ_float[islab] == aHit->getPosition()[2]) {
		i_slab = islab;
		break;
              }
	    }
	    streamlog_out(DEBUG)<<"Closest slab: "<<i_slab<<". Distance: "<<hitZtoslab.at(i_slab)<<endl;
	    	    
	    // End insertion Jesus             

	    if (_ConversionGeV2MIP) {
	      float thisenergy = aHit->getEnergy() / _GeV2MIPFactors_float[i_slab];
	      hit_energy_og.push_back(thisenergy);
	      sum_energy_og += thisenergy;
	      if(thisenergy > _MIPThreshold){
		hit_energy.push_back(thisenergy);
		sum_energy += thisenergy;
		hit_x.push_back(aHit->getPosition()[0]);
		hit_y.push_back(aHit->getPosition()[1]);
		hit_z.push_back(aHit->getPosition()[2]);
		hit_slab.push_back(i_slab);
	
		if (aHit->getPosition()[0] > 0) gap_hit_x = aHit->getPosition()[0] + _HalfCenterGap;
		else gap_hit_x = aHit->getPosition()[0] - _HalfCenterGap;
		if (aHit->getPosition()[1] > 0) gap_hit_y = aHit->getPosition()[1] + _HalfCenterGap;
		else gap_hit_y = aHit->getPosition()[1] - _HalfCenterGap;

		for (int icell = 0; icell < 1024; icell++) {
		  float in_x = fabs(_maps[_SlabMapIndices[i_slab]][icell][0] - gap_hit_x);
		  float in_y = fabs(_maps[_SlabMapIndices[i_slab]][icell][1] - gap_hit_y);
		  if (in_x < 0.1 && in_y < 0.1) {
		    hit_chip.push_back(icell/64);
		    hit_chan.push_back(icell%64);
		  }
		}
	      }	      
	    }
	    else {
	      float thisenergy = aHit->getEnergy();
	      hit_energy_og.push_back(thisenergy);
              hit_energy.push_back(thisenergy);
              sum_energy += thisenergy;
	      sum_energy_og += thisenergy;
	      hit_x.push_back(aHit->getPosition()[0]);
	      hit_y.push_back(aHit->getPosition()[1]);
	      hit_z.push_back(aHit->getPosition()[2]);
	      hit_slab.push_back(i_slab);
	      
	      if (aHit->getPosition()[0] > 0) gap_hit_x = aHit->getPosition()[0] + _HalfCenterGap;
	      else gap_hit_x = aHit->getPosition()[0] - _HalfCenterGap;
	      if (aHit->getPosition()[1] > 0) gap_hit_y = aHit->getPosition()[1] + _HalfCenterGap;
	      else gap_hit_y = aHit->getPosition()[1] - _HalfCenterGap;

	      for (int icell = 0; icell < 1024; icell++) {
		float in_x = fabs(_maps[_SlabMapIndices[i_slab]][icell][0] - gap_hit_x);
		float in_y = fabs(_maps[_SlabMapIndices[i_slab]][icell][1] - gap_hit_y);
		if (in_x < 0.1 && in_y < 0.1) {
		  hit_chip.push_back(icell/64);
		  hit_chan.push_back(icell%64);
		}
	      }
	    }

            // ** Note ** //
            // I didn't find a straightforward way to fill hit_slab,
            // without providing the list of z layer pos (_FixedPosZ) as
            // input in the steering, which is not convenient because
            // per conf/geometry, z points can be different...
            // ** End Note ** //
            
            // Shorter implementation if perfect equality between hit_z and _FixedPosZ
            // auto slab_index = std::find(_FixedPosZ.begin(), _FixedPosZ.end(), aHit->getPosition()[2]);
            // if (slab_index != _FixedPosZ.end()) hit_slab.push_back(std::distance(_FixedPosZ.begin(), slab_index));

	    // ** Note Jesus ** // 
	    // I fixed most of this issue with the insertion above
	    // Also added debugging tools to know what is going wrong in case the slabs positions is not correct
	    // or if MIP values or Si Thicknesses values are missing
	
	}//end loop over SimCalorimeterHits
     
          vector<int>::iterator it1;
          // nhit_slab
          slabs_hit = hit_slab;
          sort(slabs_hit.begin(), slabs_hit.end());
          it1 = unique(slabs_hit.begin(), slabs_hit.end());
          slabs_hit.resize(std::distance(slabs_hit.begin(), it1));
          nhit_slab = slabs_hit.size();
	  streamlog_out(DEBUG)<<"slabs_hit.size(): "<<slabs_hit.size()<<endl;

          // nhit_chan
          vector<int>::iterator it2;
          chans_hit = hit_chan;
          sort(chans_hit.begin(), chans_hit.end());
          it2 = unique(chans_hit.begin(), chans_hit.end());
          chans_hit.resize(std::distance(chans_hit.begin(), it2));
          nhit_chan = chans_hit.size();
	  streamlog_out(DEBUG)<<"chans_hit.size(): "<<chans_hit.size()<<endl;
	  
          // nhit_chip
          vector<int>::iterator it3;
          chips_hit = hit_chip;
          sort(chips_hit.begin(), chips_hit.end());
          it3 = unique(chips_hit.begin(), chips_hit.end());
          chips_hit.resize(distance(chips_hit.begin(), it3));
          nhit_chip = chips_hit.size();
	  streamlog_out(DEBUG)<<"chips_hit.size(): "<<chips_hit.size()<<endl;
	
        }//end if col
      }//end if find col names

          
    }//end for loop

    streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber()
                         << "   in run:  " << evt->getRunNumber() << std::endl ;
    if(_ConversionGeV2MIP) streamlog_out(DEBUG) <<"(Inthis event) Total hits: "<<hit_energy_og.size()<<". Hits over threshold: "<<hit_energy.size()<<std::endl;
    
    _treeout->Fill();
    // Clear vectors
    hit_slab.clear();
    hit_chip.clear();
    hit_chan.clear();
    hit_sca.clear();
    hit_adc_high.clear();
    hit_adc_low.clear();
    hit_n_scas_filled.clear();
    hit_isHit.clear();
    hit_isMasked.clear();
    hit_isCommissioned.clear();
    
    hit_energy.clear();
    hit_energy_og.clear();
    hit_x.clear();
    hit_y.clear();
    hit_z.clear();
    
    _nEvt ++ ;
    
  }

  /************************************************************************************/

  void LCIO2BuildProcessor::check( LCEvent * evt ) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
  }

  /************************************************************************************/
  void LCIO2BuildProcessor::end(){

    _treeout->Write();
    _rootout->Write("", TObject::kOverwrite);
    _rootout->Close();
    
    std::cout << "LCIO2BuildProcessor::end()  " << name()
      << " processed " << _nEvt << " events in " << _nRun << " runs "
      << "FLAG"
      << std::endl ;

  }

  /************************************************************************************/
  // void LCIO2BuildProcessor::hitCast(){
  // }
  
  // SimCalorimeterHit LCIO2BuildProcessor::hitCast(SimCalorimeterHit* aHit){
  //   return dynamic_cast<SimCalorimeterHit*>aHit;
  // }

  // SimCalorimeterHit LCIO2BuildProcessor::hitCast(SimCalorimeterHit* aHit){
  //   return dynamic_cast<SimCalorimeterHit*>aHit;
  // }

  // EVENT::CalorimeterHit* LCIO2BuildProcessor::hitCast(EVENT::CalorimeterHit* aHit){
  //   return dynamic_cast<CalorimeterHit*>aHit;
  // }

  /************************************************************************************/
  void LCIO2BuildProcessor::printParameters(){
    std::cerr << "============= LCIO2Build Processor =================" <<std::endl;
    std::cerr << "Converting Simulation Hits to build ROOT file " <<std::endl;
    std::cerr << "Input Collection name : "; for(unsigned int i = 0; i < _calorimInpCollections.size(); i++) std::cerr << _calorimInpCollections.at(i) << " ";
    std::cerr << std::endl;
    //std::cerr << "Output Collection name : "; for(unsigned int i = 0; i < _calorimOutCollections.size(); i++) std::cerr << _calorimOutCollections.at(i) << " ";
    std::cerr << std::endl;
    //std::cerr << "MIP: " << _MIPvalue << " GeV" <<std::endl;
    std::cerr << "=======================================================" <<std::endl;
    return;

  }

}
