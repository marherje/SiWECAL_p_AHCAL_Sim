#include "FLATv01Processor.hh"

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
#include "UTIL/CellIDDecoder.h"
#include "UTIL/CellIDEncoder.h"
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

  FLATv01Processor aFLATv01Processor;

  FLATv01Processor::FLATv01Processor() : Processor("FLATv01Processor")
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
    
  }

  /************************************************************************************/

  void FLATv01Processor::init()
  {
    streamlog_out(DEBUG) << "init called" << std::endl ;

    printParameters();

    gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
    gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

    _nRun = 0 ;
    _nEvt = 0 ;
    _rootout = new TFile(_eConfName.c_str(), "RECREATE");
    _treeout = new TTree("SLCIOConverted", "From SLCIO");
    // _treeout->Branch("nHits", &_nHits);
    // Hit length
    //_treeout->Branch("_nSubHits", &_nSubHits);
    _treeout->Branch("energy", &_energy);
    _treeout->Branch("positionX", &_positionX);
    _treeout->Branch("positionY", &_positionY);
    _treeout->Branch("positionZ", &_positionZ);
    _treeout->Branch("cellID0", &_cellID0);
    _treeout->Branch("cellID1", &_cellID1);
    //_treeout->Branch("_nMCParticles", &_nMCParticles);
    _treeout->Branch("nMCContributions", &_nMCContributions);
    _treeout->Branch("colType", &_colType);
    // Max subhit (hit length)
    _treeout->Branch("maxE_subhit_energy", &_maxE_subhit_energy);
    _treeout->Branch("maxE_subhit_time", &_maxE_subhit_time);
    _treeout->Branch("maxE_subhit_positionX", &_maxE_subhit_positionX);
    _treeout->Branch("maxE_subhit_positionY", &_maxE_subhit_positionY);
    _treeout->Branch("maxE_subhit_positionZ", &_maxE_subhit_positionZ);
    _treeout->Branch("maxE_subhit_PDG", &_maxE_subhit_PDG);
    _treeout->Branch("maxE_particle_pX",  &_maxE_particle_pX);
    _treeout->Branch("maxE_particle_pY",  &_maxE_particle_pY);
    _treeout->Branch("maxE_particle_pZ",  &_maxE_particle_pZ);
    _treeout->Branch("maxE_particle_energy",  &_maxE_particle_energy);
    _treeout->Branch("maxE_particle_PDG", &_maxE_particle_PDG);
    _treeout->Branch("maxE_particle_index", &_maxE_particle_index);

    _colDict["ECalBarrelSiHitsEven"]     = 1;
    _colDict["ECalBarrelSiHitsOdd"]      = 2;
    _colDict["ECalEndcapSiHitsEven"]     = 3;
    _colDict["ECalEndcapSiHitsOdd"]      = 4;
    _colDict["HcalBarrelRegCollection"]  = 5;
    _colDict["HcalEndcapRingCollection"] = 6;
    _colDict["HcalEndcapsCollection"]    = 7;

  }


  /************************************************************************************/
  void FLATv01Processor::processRunHeader( LCRunHeader* run)
  {
    _nRun++ ;
  }

  void FLATv01Processor::processEvent( LCEvent * evt )
  {

    int evtNumber = evt->getEventNumber();
    if ((evtNumber % 1000) == 0)
    streamlog_out(DEBUG) << " \n ---------> Event: "<<evtNumber<<"!Collections!! <-------------\n" << std::endl;

    //Get the list of collections in the event
    const std::vector<std::string> *cnames = evt->getCollectionNames();
    
    // Here add MC particle truth loop
    LCCollection *MCParticleCollection = evt->getCollection("MCParticle");
    int noParticles = MCParticleCollection->getNumberOfElements();
    for (int iPart = 0; iPart < noParticles; iPart++) {
      MCParticle *aParticle = dynamic_cast<MCParticle*>(MCParticleCollection->getElementAt(iPart));
      _particles.push_back(aParticle);
    }
    
    // std::cout << "collection size " << _calorimInpCollections.size() << "\n";
    //_nHits = 0;
    std::vector<float> this_energyCont;
    int this_nMC;
    for(unsigned int icol = 0; icol < _calorimInpCollections.size(); icol++)
    {
      if (icol == 2) continue;
      // std::cout << "Entering icol loop, icol = " << icol << "\n";
      if( std::find(cnames->begin(), cnames->end(), _calorimInpCollections.at(icol)) != cnames->end() )
      {
        string _inputColName = _calorimInpCollections.at(icol);
        

        //std::cout << "Entering icol if = " << _inputColName << "\n";

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
          // _nHits += noHits;
          //_nHits = inputCalorimCollection->getNumberOfElements();
          
          // Assign col type variable
          int this_colType;
          // TODO: should be assigned before and pushed back here
          if ( _colDict.find(_inputColName) == _colDict.end() ) {
            // _colType.push_back(0);
            this_colType = 0;
          } else {
            // _colType.push_back(_colDict[_inputColName]);
            this_colType = _colDict[_inputColName];
          }

          // std::cout << "Before nohits loop\n";
          for (int i = 0; i < noHits; i++)
          {
            // if (i==2){
            //   std::cout << "\tBreak at 2 hits\n";
            //   break;
            // }
            SimCalorimeterHit *aSimCalorimHit = dynamic_cast<SimCalorimeterHit*>(inputCalorimCollection->getElementAt(i));

            //_nSubHits.push_back(aSimCalorimHit->getNMCContributions());
            _energy.push_back(aSimCalorimHit->getEnergy());
            _positionX.push_back(aSimCalorimHit->getPosition()[0]);
            _positionY.push_back(aSimCalorimHit->getPosition()[1]);
            _positionZ.push_back(aSimCalorimHit->getPosition()[2]);
            _cellID0.push_back(aSimCalorimHit->getCellID0());
            _cellID1.push_back(aSimCalorimHit->getCellID1());
            //_nMCParticles.push_back(aSimCalorimHit->getNMCParticles());
            _nMCContributions.push_back(aSimCalorimHit->getNMCContributions());
            this_nMC = aSimCalorimHit->getNMCContributions();
            _colType.push_back(this_colType);
            
            //std::cout << "\tHit " << i << " values pushed back\n";

            //std::cout << "\tthis_energyCont before fill size " << this_energyCont.size() << "\n";

            //this_energyCont.resize(_nMCContributions[i], 0.0);
            this_energyCont.resize(this_nMC, 0.0);
            // for (int it = 0; it < _nMCContributions[i]; it++)
            for (int it = 0; it < this_nMC; it++)
            {
              //std::cout << "\t\tSubhit cont " << it << " of " << _nMCContributions[i] << ", " << aSimCalorimHit->getEnergyCont(it) << " GeV\n";
              // this_energyCont.push_back(aSimCalorimHit->getEnergyCont(it));
              this_energyCont[it] = aSimCalorimHit->getEnergyCont(it);
            }
            // std::cout << "Between print and push\n";
            // for (int it = 0; it < _nMCContributions[i]; it++)
            // {
            //   // this_energyCont.push_back(aSimCalorimHit->getEnergyCont(it));
            //   this_energyCont[it] = aSimCalorimHit->getEnergyCont(it);
            // }
            
            // std::cout << "\tthis_energyCont after fill size " << this_energyCont.size() << "\n";
            
            int argMaxE = std::distance(this_energyCont.begin(), std::max_element(this_energyCont.begin(), this_energyCont.end()));

            _maxE_subhit_energy.push_back(this_energyCont[argMaxE]);
            _maxE_subhit_time.push_back(aSimCalorimHit->getTimeCont(argMaxE));
            _maxE_subhit_positionX.push_back(aSimCalorimHit->getStepPosition(argMaxE)[0]);
            _maxE_subhit_positionY.push_back(aSimCalorimHit->getStepPosition(argMaxE)[1]);
            _maxE_subhit_positionZ.push_back(aSimCalorimHit->getStepPosition(argMaxE)[2]);
            _maxE_subhit_PDG.push_back(aSimCalorimHit->getPDGCont(argMaxE));

            auto maxPart = aSimCalorimHit->getParticleCont(argMaxE);
            _maxE_particle_pX.push_back(maxPart->getMomentum()[0]);
            _maxE_particle_pY.push_back(maxPart->getMomentum()[1]);
            _maxE_particle_pZ.push_back(maxPart->getMomentum()[2]);
            _maxE_particle_energy.push_back(maxPart->getEnergy());
            _maxE_particle_PDG.push_back(maxPart->getPDG());
            
			// Find MCParticle index in event
            auto it = std::find(_particles.begin(), _particles.end(), maxPart);
            if (it != _particles.end() ){
              _maxE_particle_index.push_back(it - _particles.begin());
            }
            else _maxE_particle_index.push_back(-1);

            //this_energyCont.clear();

            //std::cout << "Break SimCalorimeterHits";
            // break;
          }//end loop over SimCalorimeterHits
          //std::cout << "End of hit loop"<< endl;
        // Used to be at this level: _treeout->Fill();
        }//end if col
        
        
      }//end if find col names
      //std::cout << "Break colnames";

          
      // std::cout << "Ended first loop, icol = " << icol << "\n";
      // std::cout << "_maxE_particle_pX size = " << _maxE_particle_pX.size() << "\n";
      // std::cout << "_maxE_particle_PDG size = " << _maxE_particle_PDG.size() << "\n";
      // std::cout << "_maxE_subhit_energy size = " << _maxE_subhit_energy.size() << "\n";
      // if (icol == 3) break;
    }//end for loop

  // std::cout << "After loop, before fill";
  _treeout->Fill();
  //_nSubHits.clear();
  // Clear vectors
  _energy.clear();
  _positionX.clear();
  _positionY.clear();
  _positionZ.clear();
  _cellID0.clear();
  _cellID1.clear();
  //_nMCParticles.clear();
  _nMCContributions.clear();
  _colType.clear();
  //
  _maxE_subhit_time.clear();
  _maxE_subhit_energy.clear();
  _maxE_subhit_positionX.clear();
  _maxE_subhit_positionY.clear();
  _maxE_subhit_positionZ.clear();
  _maxE_subhit_PDG.clear();
  _maxE_particle_pX.clear();
  _maxE_particle_pY.clear();
  _maxE_particle_pZ.clear();
  _maxE_particle_energy.clear();
  _maxE_particle_PDG.clear();
  _maxE_particle_index.clear();
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber()
  << "   in run:  " << evt->getRunNumber() << std::endl ;


  _nEvt ++ ;

  // if (_nEvt % 1000 == 0) {
  //   _treeout->Write("", TObject::kOverwrite);
  //   _rootout->Write("", TObject::kOverwrite);
  // }

  _particles.clear();
  }

  /************************************************************************************/

  void FLATv01Processor::check( LCEvent * evt ) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
  }

  /************************************************************************************/
  void FLATv01Processor::end(){

    _treeout->Write();
    _rootout->Write("", TObject::kOverwrite);
    _rootout->Close();
    
    std::cout << "FLATv01Processor::end()  " << name()
      << " processed " << _nEvt << " events in " << _nRun << " runs "
      << "FLAG"
      << std::endl ;

  }


  /************************************************************************************/
  void FLATv01Processor::printParameters(){
    std::cerr << "============= FLATv01 Processor =================" <<std::endl;
    //std::cerr << "Converting Simulation Hits from GeV to MIP" <<std::endl;
    std::cerr << "Converting Simulation Hits to flat ROOT file (v01)" <<std::endl;
    std::cerr << "Input Collection name : "; for(unsigned int i = 0; i < _calorimInpCollections.size(); i++) std::cerr << _calorimInpCollections.at(i) << " ";
    std::cerr << std::endl;
    //std::cerr << "Output Collection name : "; for(unsigned int i = 0; i < _calorimOutCollections.size(); i++) std::cerr << _calorimOutCollections.at(i) << " ";
    std::cerr << std::endl;
    //std::cerr << "MIP: " << _MIPvalue << " GeV" <<std::endl;
    std::cerr << "=======================================================" <<std::endl;
    return;

  }


}
