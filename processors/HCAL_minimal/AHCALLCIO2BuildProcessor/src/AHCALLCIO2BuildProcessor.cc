#include "AHCALLCIO2BuildProcessor.hh"

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

  AHCALLCIO2BuildProcessor aAHCALLCIO2BuildProcessor;

  AHCALLCIO2BuildProcessor::AHCALLCIO2BuildProcessor() : Processor("AHCALLCIO2BuildProcessor")
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

    vector<string> FixedPosZExample;
    registerProcessorParameter("FixedPosZ",
                               "vector of z slabs positions",
                               _FixedPosZ,
                               FixedPosZExample);
        
    registerProcessorParameter("NSlabs",
                               "Number of slabs",
                               _NSlabs,
                               int(9999));

    registerProcessorParameter("FirstSlabsPosZ",
                               "Position of first slab",
                               _FirstSlabPosZ,
                               float(9999));

    registerProcessorParameter("SlabSpacing",
                               "Separation between slabs",
                               _SlabSpacing,
                               float(9999));
    
    registerProcessorParameter("SlabThickness",
                               "Scintillator thickness per slab",
                               _SlabThickness,
                               float(9999));

    registerProcessorParameter("ConversionGeV2MIP",
                               "Do conversion to MIP",
                               _ConversionGeV2MIP,
                               false);

    registerProcessorParameter("GeV2MIPFactor",
                               "GeV2MIP conversion",
			       _GeV2MIPFactor,
                               float(9999));

    registerProcessorParameter("MIPThreshold",
                               "Minimum energy (in MIPs) to count a hit",
                               _MIPThreshold,
                               float(0.));
            
  }

  /************************************************************************************/

  void AHCALLCIO2BuildProcessor::init()
  {
    streamlog_out(DEBUG) << "init called" << std::endl ;

    printParameters();

    _nRun = 0 ;
    _nEvt = 0 ;
    _rootout = new TFile(_outFileName.c_str(), "RECREATE");
    _treeout = new TTree("hcal", "From SLCIO");
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
      for (int islab = 0; islab < int(_FixedPosZ.size()); islab++){
        _FixedPosZ_float.push_back(stof(_FixedPosZ[islab]));
      }
      if(_SlabThickness == 9999.){ streamlog_out(ERROR)<<"Missing slab thickness"<<std::endl;}
    }
    else{ 
      if(_NSlabs == 9999){ streamlog_out(ERROR)<<"Missing number of slabs"<<std::endl;}
      else {
	if(_FirstSlabPosZ == 9999.){ streamlog_out(ERROR)<<"Missing first slab position"<<std::endl;}
	else {
	  if(_SlabSpacing == 9999.){ streamlog_out(ERROR)<<"Missing slab separation"<<std::endl;}
	  else {
	    if(_SlabThickness == 9999.){ streamlog_out(ERROR)<<"Missing slab thickness"<<std::endl;}
	    else{
	      for (int islab = 0; islab < _NSlabs; islab++){
		_FixedPosZ_float.push_back(_FirstSlabPosZ + _SlabSpacing*islab);
	      }
	    }
	  }
	}
      }
    }
   

    if(_ConversionGeV2MIP==true){
      if(_GeV2MIPFactor==9999) streamlog_out(ERROR)<< "Missing MIP info "<<std::endl;
    }
    
    streamlog_out(DEBUG)<<"PARAMETER READING"<<std::endl;
    streamlog_out(DEBUG)<<"GeV2MIP: "<<_GeV2MIPFactor<<std::endl;
    streamlog_out(DEBUG)<<"Threshold energy (MIPS): "<<_MIPThreshold<<std::endl;
    streamlog_out(DEBUG)<<"FixedPosZ size: "<<_FixedPosZ.size()<<std::endl;
    streamlog_out(DEBUG)<<"N Slabs: "<<_NSlabs<<", First slab pos Z: "<<_FirstSlabPosZ<<", Slab spacing: "<<_SlabSpacing<<std::endl;
    streamlog_out(DEBUG)<<"Detector Thickness: "<<_SlabThickness<<std::endl;

    _printType = true;

  }


  /************************************************************************************/
  void AHCALLCIO2BuildProcessor::processRunHeader( LCRunHeader* run)
  {
    _nRun++ ;
  }

  void AHCALLCIO2BuildProcessor::processEvent( LCEvent * evt )
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
          int i_slab = 0;
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

	    vector<float> hitZtoslab;
	    float smallestdistance=9999.;
	    float hitZ=aHit->getPosition()[2];
	    
	    if( _FixedPosZ_float.size() == 0) streamlog_out(ERROR)<<"missing layer info"<<endl;
	    	    
	    if (_ConversionGeV2MIP) {
              float thisenergy = aHit->getEnergy() / _GeV2MIPFactor;
              hit_energy_og.push_back(thisenergy);
              sum_energy_og += thisenergy;
              if(thisenergy > _MIPThreshold){
                hit_energy.push_back(thisenergy);
                sum_energy += thisenergy;
                hit_x.push_back(aHit->getPosition()[0]);
                hit_y.push_back(aHit->getPosition()[1]);
                hit_z.push_back(aHit->getPosition()[2]);
              
		for (int islab = 0; islab < _FixedPosZ_float.size(); islab++){
		  hitZtoslab.push_back(abs(_FixedPosZ_float[islab]-hitZ));
		}
		for (int islab = 0; islab < _FixedPosZ_float.size(); islab++){
		  streamlog_out(DEBUG)<<"layer "<<islab<<" distance: "<<hitZtoslab.at(islab)<<endl;
		  if(hitZtoslab.at(islab) < _SlabThickness) {
		    i_slab = islab;
		    hit_slab.push_back(i_slab);
		    break;
		  }
		}
		streamlog_out(DEBUG)<<"Closest slab: "<<i_slab<<". Distance: "<<hitZtoslab.at(i_slab)<<endl;
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
            
	      for (int islab = 0; islab < _FixedPosZ_float.size(); islab++){
		hitZtoslab.push_back(abs(_FixedPosZ_float[islab]-hitZ));
	      }
	      for (int islab = 0; islab < _FixedPosZ_float.size(); islab++){
		streamlog_out(DEBUG)<<"layer "<<islab<<" distance: "<<hitZtoslab.at(islab)<<endl;
		if(hitZtoslab.at(islab) < _SlabThickness) {
		  i_slab = islab;
		  hit_slab.push_back(i_slab);
		  break;
		}
	      }
	      streamlog_out(DEBUG)<<"Closest slab: "<<i_slab<<". Distance: "<<hitZtoslab.at(i_slab)<<endl;
	    }
	    
	    
	  }//end loop over SimCalorimeterHits

        }//end if col
      }//end if find col names

          
    }//end for loop

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
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber()
  << "   in run:  " << evt->getRunNumber() << std::endl ;


  _nEvt ++ ;

  }

  /************************************************************************************/

  void AHCALLCIO2BuildProcessor::check( LCEvent * evt ) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
  }

  /************************************************************************************/
  void AHCALLCIO2BuildProcessor::end(){

    _treeout->Write();
    _rootout->Write("", TObject::kOverwrite);
    _rootout->Close();
    
    std::cout << "AHCALLCIO2BuildProcessor::end()  " << name()
      << " processed " << _nEvt << " events in " << _nRun << " runs "
      << "FLAG"
      << std::endl ;

  }

  /************************************************************************************/
  // void AHCALLCIO2BuildProcessor::hitCast(){
  // }
  
  // SimCalorimeterHit AHCALLCIO2BuildProcessor::hitCast(SimCalorimeterHit* aHit){
  //   return dynamic_cast<SimCalorimeterHit*>aHit;
  // }

  // SimCalorimeterHit AHCALLCIO2BuildProcessor::hitCast(SimCalorimeterHit* aHit){
  //   return dynamic_cast<SimCalorimeterHit*>aHit;
  // }

  // EVENT::CalorimeterHit* AHCALLCIO2BuildProcessor::hitCast(EVENT::CalorimeterHit* aHit){
  //   return dynamic_cast<CalorimeterHit*>aHit;
  // }

  /************************************************************************************/
  void AHCALLCIO2BuildProcessor::printParameters(){
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
