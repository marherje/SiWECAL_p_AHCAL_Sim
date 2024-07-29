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
                               "vector of z layer positions",
                               _FixedPosZ,
                               FixedPosZExample);
        
    registerProcessorParameter("NLayers",
                               "Number of layers",
                               _NLayers,
                               int(9999));

    registerProcessorParameter("FirstLayerPosZ",
                               "Position of first layer",
                               _FirstLayerPosZ,
                               float(9999));

    registerProcessorParameter("LayerSpacing",
                               "Separation between layer",
                               _LayerSpacing,
                               float(9999));
    
    registerProcessorParameter("LayerThickness",
                               "Thickness of one layer",
                               _LayerThickness,
                               float(9999));

    registerProcessorParameter("ConversionGeV2MIP",
                               "Do conversion to MIP",
                               _ConversionGeV2MIP,
                               false);

    registerProcessorParameter("GeV2MIPFactor_float",
                               "GeV2MIP conversion",
			       _GeV2MIPFactor_float,
                               float(9999));
            
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
    _treeout->Branch("nhit_layer", &nhit_layer);
    _treeout->Branch("nhit_chip", &nhit_chip);
    _treeout->Branch("nhit_chan", &nhit_chan);
    _treeout->Branch("nhit_len", &nhit_len);

    // _treeout->Branch("sum_hg", &sum_hg);
    _treeout->Branch("sum_energy", &sum_energy);
    _treeout->Branch("sum_energy_lg", &sum_energy_lg);

    _treeout->Branch("hit_layer", &hit_layer);
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
    _treeout->Branch("hit_energy_lg", &hit_energy_lg);
    _treeout->Branch("hit_x", &hit_x);
    _treeout->Branch("hit_y", &hit_y);
    _treeout->Branch("hit_z", &hit_z);


    if(_FixedPosZ.size() != 0) {
      for (int ilayer = 0; ilayer < int(_FixedPosZ.size()); ilayer++){
        _FixedPosZ_float.push_back(stof(_FixedPosZ[ilayer]));
      }
    }
    else{ 
      if(_NLayers == 9999){ streamlog_out(ERROR)<<"Missing number of layers"<<std::endl;}
      else {
	if(_FirstLayerPosZ == 9999.){ streamlog_out(ERROR)<<"Missing first layer position"<<std::endl;}
	else {
	  if(_LayerSpacing == 9999.){ streamlog_out(ERROR)<<"Missing layer separation"<<std::endl;}
	  else {
	    if(_LayerThickness == 9999.){ streamlog_out(ERROR)<<"Missing layer thickness"<<std::endl;}
	    else{
	      for (int ilayer = 0; ilayer < _NLayers; ilayer++){
		_FixedPosZ_float.push_back(_FirstLayerPosZ + _LayerSpacing*ilayer);
	      }
	    }
	  }
	}
      }
    }
   

    if(_ConversionGeV2MIP==true){
      if(_GeV2MIPFactor_float==9999) streamlog_out(ERROR)<< "Missing MIP info "<<std::endl;
    }
    
    streamlog_out(DEBUG)<<"GeV2MIP value: "<<_GeV2MIPFactor_float<<std::endl;
    streamlog_out(DEBUG)<<"FixedPosZ size: "<<_FixedPosZ.size()<<std::endl;
    streamlog_out(DEBUG)<<"N Layers: "<<_NLayers<<", First layer pos Z: "<<_FirstLayerPosZ<<", Layer spacing: "<<_LayerSpacing<<std::endl;

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

          sum_energy_lg = 0.;
          sum_energy = 0.;

          vector<int> layers_hit, chans_hit, chips_hit;
          int i_layer;
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

            hit_energy_lg.push_back(aHit->getEnergy());
            hit_x.push_back(aHit->getPosition()[0]);
            hit_y.push_back(aHit->getPosition()[1]);
            hit_z.push_back(aHit->getPosition()[2]);
    
	    //streamlog_out(DEBUG)<<"CellID0: "<<aHit->getCellID0()<<", CellID1: "<<aHit->getCellID1()<<endl;

	    // Insertion Jesus                                                                  
	    vector<float> hitZtolayer;
	    float smallestdistance=9999.;
	    float hitZ=aHit->getPosition()[2];
	    int NLayers=0;

            if(_NLayers != 9999.) { NLayers=_NLayers; }
            else { streamlog_out(ERROR)<<"Missing number of layers"<<std::endl; }

	    for (int ilayer = 0; ilayer < NLayers; ilayer++){
	      hitZtolayer.push_back(abs(_FixedPosZ_float[ilayer]-hitZ));
	      if(ilayer==0) {
                smallestdistance=hitZtolayer.at(0);
                i_layer=0;
              }
              else if(hitZtolayer.at(ilayer) == 0.) {
                smallestdistance=hitZtolayer.at(ilayer);
                i_layer=ilayer;
                break;
              }
              else if(hitZtolayer.at(ilayer) < smallestdistance) {
                smallestdistance=hitZtolayer.at(ilayer);
                i_layer=ilayer;
              }
            }
            streamlog_out(DEBUG)<<"Closest layer: "<<i_layer<<". Distance: "<<hitZtolayer.at(i_layer)<<endl;
	    
	    for (int ilayer = 0; ilayer < NLayers; ilayer++){
	      float thislayerpos = _FirstLayerPosZ + ilayer*_LayerSpacing;
	      if ( (thislayerpos > (aHit->getPosition()[2] - _LayerThickness) ) &&
		(thislayerpos < (aHit->getPosition()[2] + _LayerThickness) ) ) {
		hit_layer.push_back(ilayer);
	      }
	      else if (thislayerpos == aHit->getPosition()[2]) {
		hit_layer.push_back(ilayer);
	      }
	    }
	    
	    // End insertion Jesus             

	    if (_ConversionGeV2MIP) {
              hit_energy.push_back(aHit->getEnergy() / _GeV2MIPFactor_float);
	      sum_energy += aHit->getEnergy() / _GeV2MIPFactor_float;
            }
            else {
              hit_energy.push_back(aHit->getEnergy());
              sum_energy += aHit->getEnergy();
            }
          }//end loop over SimCalorimeterHits

          sum_energy_lg = sum_energy;

        }//end if col
      }//end if find col names

          
    }//end for loop

  _treeout->Fill();
  // Clear vectors
  hit_layer.clear();
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
  hit_energy_lg.clear();
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
