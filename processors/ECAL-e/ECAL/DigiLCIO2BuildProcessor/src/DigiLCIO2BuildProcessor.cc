#include"DigiLCIO2BuildProcessor.hh"

#include <iostream>
#include <cstdlib>
#include <map>

// ---- LCIO Headers
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include "IMPL/CalorimeterHitImpl.h"
#include "IMPL/SimCalorimeterHitImpl.h"
#include <IMPL/LCFlagImpl.h>
#include "UTIL/LCTypedVector.h"
// #include "UTIL/CellIDDecoder.h"
// #include "UTIL/CellIDEncoder.h"
#include "UTIL/BitField64.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
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

  DigiLCIO2BuildProcessor aDigiLCIO2BuildProcessor;

  DigiLCIO2BuildProcessor::DigiLCIO2BuildProcessor() : Processor("DigiLCIO2BuildProcessor")
  {

    vector<string> calorimInpCollectionsExample;
    registerProcessorParameter("Input_Collections",
                               "vector of collections for input",
                               _calorimInpCollections,
                               calorimInpCollectionsExample);

    string eConfNameExample;
    registerProcessorParameter("Energy_Conf_Name",
                               "Name to identify output file",
                               _eConfName,
                               eConfNameExample);
    
    vector<string> siThicknessesExample;
    registerProcessorParameter("SiThicknesses",
                               "vector of Silicon thicknesses per layer",
                               _siThicknesses,
                               siThicknessesExample);
        
    // Two methods
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
    
    // registerProcessorParameter("hitType",
    //                            "Hit type (SimCalorimeterHit or CalorimeterHit)",
    //                            _hitType,
    //                            "SimCalorimeterHit");
    
  }

  /************************************************************************************/

  void DigiLCIO2BuildProcessor::init()
  {
    streamlog_out(DEBUG) << "init called" << std::endl ;

    printParameters();

    _nRun = 0 ;
    _nEvt = 0 ;
    _rootout = new TFile(_eConfName.c_str(), "RECREATE");
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
    _treeout->Branch("sum_energy_lg", &sum_energy_lg);

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
    _treeout->Branch("hit_energy_lg", &hit_energy_lg);
    _treeout->Branch("hit_x", &hit_x);
    _treeout->Branch("hit_y", &hit_y);
    _treeout->Branch("hit_z", &hit_z);
    // _treeout->Branch("cellID0", &_cellID0);
    // _treeout->Branch("cellID1", &_cellID1);
    
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
	  _siThicknesses_float.push_back(0.001*stof(_siThicknesses[ilayer])); // Note the conversion from micra to mm
	}
      }
      else {
	try{
	  if(_siThicknesses.size() == 1) _siThicknesses_float.push_back(0.001*stof(_siThicknesses[0])); // Note the conversion from micra to mm
	  else throw 404;
	}
	catch(...){streamlog_out(ERROR)<< "Missing Si Thicknesses info "<<std::endl;
	}
      }
    }

    _printType = true;
  }
  
  



  /************************************************************************************/
  void DigiLCIO2BuildProcessor::processRunHeader( LCRunHeader* run)
  {
    _nRun++ ;
  }

  void DigiLCIO2BuildProcessor::processEvent( LCEvent * evt )
  {
    
    int evtNumber = evt->getEventNumber();
    if ((evtNumber % 1000) == 0)
    streamlog_out(DEBUG) << " \n ---------> Event: "<<evtNumber<<"!Collections!! <-------------\n" << std::endl;

    //Get the list of collections in the event
    const std::vector<std::string> *cnames = evt->getCollectionNames();
    
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
          nhit_chip = _FixedPosZ.size() * 16;
          nhit_chan = _FixedPosZ.size() * 16 * 64;
          nhit_len = noHits;
          // nhit_slab = 0;

          // sum_hg = -1;
          // sum_energy_lg = -1;

          sum_energy = 0.;

          for (int i = 0; i < noHits; i++)
          {
            CalorimeterHit *aHit = dynamic_cast<CalorimeterHit*>(inputCalorimCollection->getElementAt(i));
            //auto *aHit = hitCast(inputCalorimCollection->getElementAt(i));
            // hitCast(*aHit);

            hit_chip.push_back(-1);
            hit_chan.push_back(-1);
            hit_sca.push_back(1);
            hit_adc_high.push_back(-1);
            hit_adc_low.push_back(-1);
            hit_n_scas_filled.push_back(1);
            hit_isHit.push_back(1);
            hit_isMasked.push_back(0);
            hit_isCommissioned.push_back(1);

            hit_energy.push_back(aHit->getEnergy());
            sum_energy += aHit->getEnergy();
            hit_energy_lg.push_back(aHit->getEnergy());
            hit_x.push_back(aHit->getPosition()[0]);
            hit_y.push_back(aHit->getPosition()[1]);
            hit_z.push_back(aHit->getPosition()[2]);
            
            // Need Si Thickness tolerance to find layer
            for (int i_slab = 0; i_slab < int(_FixedPosZ.size()); i_slab++){
              if (_FixedPosZ_float[i_slab] > (aHit->getPosition()[2] - _siThicknesses_float[i_slab]) &&
                  _FixedPosZ_float[i_slab] < (aHit->getPosition()[2] + _siThicknesses_float[i_slab]) ) {
                hit_slab.push_back(i_slab);
              }
	      else if (_FixedPosZ_float[i_slab] == aHit->getPosition()[2]) {
		hit_slab.push_back(i_slab);
	      }
	    }
	    
            // Shorter implementation if perfect equality between hit_z and _FixedPosZ
            // auto slab_index = std::find(_FixedPosZ.begin(), _FixedPosZ.end(), aHit->getPosition()[2]);
            // if (slab_index != _FixedPosZ.end()) hit_slab.push_back(std::distance(_FixedPosZ.begin(), slab_index));

          }//end loop over SimCalorimeterHits

          sum_energy_lg = sum_energy;

          std::vector<float> slabs_hit;
          slabs_hit = hit_z;
          
          std::sort(slabs_hit.begin(), slabs_hit.end());
          
          std::vector<float>::iterator it;
          it = std::unique(slabs_hit.begin(), slabs_hit.end());
          slabs_hit.resize(std::distance(slabs_hit.begin(),it));
          nhit_slab = slabs_hit.size();

          //std::cout << "End of hit loop"<< endl;
        }//end if col
        
        
      }//end if find col names
      //std::cout << "Break colnames";

          
    }//end for loop

  // std::cout << "After loop, before fill";
  _treeout->Fill();
  //_nSubHits.clear();
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
  hit_energy_lg.clear();
  hit_x.clear();
  hit_y.clear();
  hit_z.clear();
  //_cellID0.clear();
  // _cellID1.clear();
  //_nMCParticles.clear();
  // _nMCContributions.clear();
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber()
  << "   in run:  " << evt->getRunNumber() << std::endl ;


  _nEvt ++ ;

  }

  /************************************************************************************/

  void DigiLCIO2BuildProcessor::check( LCEvent * evt ) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
  }

  /************************************************************************************/
  void DigiLCIO2BuildProcessor::end(){

    _treeout->Write();
    _rootout->Write("", TObject::kOverwrite);
    _rootout->Close();
    
    std::cout << "DigiLCIO2BuildProcessor::end()  " << name()
      << " processed " << _nEvt << " events in " << _nRun << " runs "
      << "FLAG"
      << std::endl ;

  }

  /************************************************************************************/
  // void DigiLCIO2BuildProcessor::hitCast(){
  // }
  
  // SimCalorimeterHit DigiLCIO2BuildProcessor::hitCast(SimCalorimeterHit* aHit){
  //   return dynamic_cast<SimCalorimeterHit*>aHit;
  // }

  // SimCalorimeterHit DigiLCIO2BuildProcessor::hitCast(SimCalorimeterHit* aHit){
  //   return dynamic_cast<SimCalorimeterHit*>aHit;
  // }

  // EVENT::CalorimeterHit* DigiLCIO2BuildProcessor::hitCast(EVENT::CalorimeterHit* aHit){
  //   return dynamic_cast<CalorimeterHit*>aHit;
  // }

  /************************************************************************************/
  void DigiLCIO2BuildProcessor::printParameters(){
    std::cerr << "============= LCIO2Build Processor =================" <<std::endl;
    //std::cerr << "Converting Simulation Hits from GeV to MIP" <<std::endl;
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
