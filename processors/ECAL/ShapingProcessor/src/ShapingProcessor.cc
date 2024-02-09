#include "ShapingProcessor.hh"
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>


#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include "IMPL/CalorimeterHitImpl.h"
#include "IMPL/SimCalorimeterHitImpl.h"
#include <IMPL/LCFlagImpl.h>
#include "EVENT/CalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "UTIL/LCTypedVector.h"
#include "UTIL/CellIDDecoder.h"
#include "UTIL/CellIDEncoder.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include "marlin/Exceptions.h"

// ----- ROOT --------
#include "TMath.h"
#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TTreeReader.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;

namespace CALICE
{

  ShapingProcessor aShapingProcessor ;


  ShapingProcessor::ShapingProcessor() : Processor("ShapingProcessor") {

    // register steering parameters: name, description, class-variable, default value

    vector<string> in_calorimInpCollections;
    registerProcessorParameter("Input_Collections",
                               "vector of collections for input",
                               _calorimInpCollections,
                               in_calorimInpCollections);

    vector<string> in_calorimOutCollections;
    registerProcessorParameter("Output_Collections",
                               "vector of collections for output",
                               _calorimOutCollections,
                               in_calorimOutCollections);
    
    // Parameters from steering

    registerProcessorParameter("ShapingProcessor_nbinsI",
                               "Number of bins of input hit histogram",
                               _nbinsI,
                               int(400));
    
    registerProcessorParameter("ShapingProcessor_nbinsF",
                               "Number of bins of Fast Shaper histogram",
                               _nbinsF,
                               int(400));

    registerProcessorParameter("ShapingProcessor_nbinsS",
                               "Number of bins of Slow Shaper histogram",
                               _nbinsS,
                               int(100));

    registerProcessorParameter("ShapingProcessor_bwI",
                               "Bin width of input hit histogram",
                               _bwI,
                               double(0.5));

    registerProcessorParameter("ShapingProcessor_bwF",
                               "Bin width of Fast Shaper histogram",
                               _bwF,
                               double(0.5));

    registerProcessorParameter("ShapingProcessor_bwS",
                               "Bin width of Slow Shaper histogram",
                               _bwS,
                               double(5));

    registerProcessorParameter("ShapingProcessor_delay",
                               "Delay to read Slow Shaper after threshold passed",
                               _delay,
                               double(160));
    
    registerProcessorParameter("ShapingProcessor_MIPThreshold",
                               "Threshold for Fast Shaper in MIP units",
                               _MIPThreshold,
                               double(0.3333));

    registerProcessorParameter("ShapingProcessor_FSNoise",
                               "Fast Shaper noise",
                               _FSNoise,
                               double(1/12));
                               //double(1));

    registerProcessorParameter("ShapingProcessor_SSNoise",
                               "Slow Shaper noise",
                               _SSNoise,
                               double(1/20));
                               //double(1/5));

    registerProcessorParameter("ShapingProcessor_useHistInput",
                               "Use histogrammed hits as input for shaper",
                               _useHistInput,
                               true);

    registerProcessorParameter("ShapingProcessor_filterNoise",
                               "Add noise to Fast and Slow Shapers",
                               _filterNoise,
                               true);

    string in_shapingAuxFilename = "aux_shaping.root";
    registerProcessorParameter("ShapingProcessor_AuxFilename",
                               "Name to identify auxiliary file",
                               _shapingAuxFilename,
                               in_shapingAuxFilename);

    
    // "Fixed" paremeters:
    // Reading MIP conversion
    vector<string> MIP2GeVFactorsExample;
    registerProcessorParameter("MIP2GeVFactors",
                               "vector of conversion factors",
                               _MIP2GeVFactors,
                               MIP2GeVFactorsExample);
    
    // Reading positions (2 methods) 
    registerProcessorParameter("deltaZ",
                               "Tolerance (mm) for assigning hits to a slab",
                               _deltaZ,
                               float(0.1));

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

    vector<string> FixedPosZExample;
    registerProcessorParameter("FixedPosZ",
                               "vector of z layer positions",
                               _FixedPosZ,
                               FixedPosZExample);
    
    // Shapers parameter
    registerProcessorParameter("tauF",
                               "Tau parameter of Fast Shaper",
                               _tauF,
                               double(30));
    
    registerProcessorParameter("tauS",
                               "Tau parameter of Slow Shaper",
                               _tauS,
                               double(180));

    registerProcessorParameter("nF",
                               "Order of CR-RC filter for Fast Shaper",
                               _nF,
                               int(2));
    
    registerProcessorParameter("nS",
                               "Order of CR-RC filter for Slow Shaper",
                               _nS,
                               int(2));
  }

  /************************************************************************************/

  void ShapingProcessor::init() {

    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _nRun = 0 ;
    _nEvt = 0 ;
    _nLeft= 0 ;
    _nTot = 0 ;

    _plotHit = false;
    _plotted = false;

    // No idea why this was here (Jesus)
    //_threshold = _MIP2GeV * _MIPThreshold;   
    //cout << "Th is mip*" << to_string(_MIPThreshold) << " = " <<  to_string(_threshold) << endl;

    streamlog_out(DEBUG)<<"MIP2GeV size: "<<_MIP2GeVFactors.size()<<std::endl;
    streamlog_out(DEBUG)<<"FixedPosZ size: "<<_FixedPosZ.size()<<std::endl;
    streamlog_out(DEBUG)<<"N Slabs: "<<_NSlabs<<", First slab pos Z: "<<_FirstSlabPosZ<<", Slab spacing: "<<_SlabSpacing<<std::endl;
    
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
      if( _MIP2GeVFactors.size() > 1 ) {
	if(_FixedPosZ.size() != _MIP2GeVFactors.size()) {streamlog_out(ERROR)<< "MIP vector size don't match the number of layers "<<std::endl;}
	else {
	  _MIP2GeVFactors_float.push_back(stof(_MIP2GeVFactors[ilayer]));
	}
      }
      else {
	try{
	  if(_MIP2GeVFactors.size() == 1) _MIP2GeVFactors_float.push_back(stof(_MIP2GeVFactors.at(0)));
	  else throw 404;
	}
	catch(...){streamlog_out(ERROR)<< "Missing MIP info "<<std::endl;
	}
      }
    }

    streamlog_out(DEBUG)<<"MIP2GeV size: "<<_MIP2GeVFactors.size()<<std::endl;
    streamlog_out(DEBUG)<<"FixedPosZ size: "<<_FixedPosZ.size()<<std::endl;
    streamlog_out(DEBUG)<<"Tolerance read (mm): "<<_deltaZ<<std::endl;
    streamlog_out(DEBUG)<<"N Slabs: "<<_NSlabs<<", First slab pos Z: "<<_FirstSlabPosZ<<", Slab spacing: "<<_SlabSpacing<<std::endl;

    _rootout = new TFile(_shapingAuxFilename.c_str(), "RECREATE");
    _treeout = new TTree("CalHits","Calorimeter hits");
    _treeout->Branch("posx", &posx, "posx/D");
    _treeout->Branch("posy", &posy, "posy/D");
    _treeout->Branch("posz", &posz, "posz/D");
    _treeout->Branch("energy_sim", &energy_sim, "energy_sim/D");
    _treeout->Branch("time_sim", &time_sim, "time_sim/D");
    _treeout->Branch("energy_sha", &energy_sha, "energy_sha/D");
    _treeout->Branch("time_sha", &time_sha, "time_sha/D");
    _treeout->Branch("t_maxFS", &t_maxFS, "t_maxFS/D");
    _treeout->Branch("t_maxSS", &t_maxSS, "t_maxSS/D");
    _treeout->Branch("maxFS", &maxFS, "maxFS/D");
    _treeout->Branch("maxSS", &maxSS, "maxSS/D");

    printParameters() ;
  }

  /************************************************************************************/
  void ShapingProcessor::processRunHeader( LCRunHeader* run) {

    _nRun++ ;
  }

  /************************************************************************************/

  void ShapingProcessor::processEvent( LCEvent * evt ) {

    int evtNumber = evt->getEventNumber();

    if(evtNumber%1000 == 0)
    streamlog_out(MESSAGE) << " \n ---------> Event: " << evtNumber <<"!!! <-------------\n"<< endl;

    //Get the list of collections in the event
    const std::vector<std::string> *cnames = evt->getCollectionNames();

    for(unsigned int icol = 0; icol < _calorimInpCollections.size(); icol++)
    {
      if( std::find(cnames->begin(), cnames->end(), _calorimInpCollections.at(icol)) != cnames->end() )
      {
        string _inputColName = _calorimInpCollections.at(icol);
        // cout << "Look at " << _inputColName << endl;

        LCCollection *inputCol = 0;
        try
        {
          inputCol = evt->getCollection(_inputColName);
        }
        catch (EVENT::DataNotAvailableException &e)
        {
	  cout<<"collection error"<<endl;
          streamlog_out(WARNING)<< "missing collection "
				<<_inputColName<<endl<<e.what()<<std::endl;
        }

        if (inputCol != 0){

          LCCollectionVec* calOutVec = new LCCollectionVec( LCIO::CALORIMETERHIT) ;
          LCFlagImpl hitFlag(calOutVec->getFlag());
          hitFlag.setBit(LCIO::RCHBIT_TIME);
          hitFlag.setBit(LCIO::CHBIT_LONG);
          calOutVec->setFlag(hitFlag.getFlag());
	  
	  //Insertion JP start
	  //LCParameters &ogParam = inputCol->parameters();                                                                                                                                                  
          //string InColEncoder=ogParam.getStringVal("LCIO::CellIDEncoding");                                                                                                     
          //LCParameters &theParam = calOutVec->parameters();                                                                                                                         
          //theParam.setValue(LCIO::CellIDEncoding, InColEncoder);
	  //Insertion JP end

          streamlog_out(DEBUG) << "Collection " << _inputColName << " found " << std::endl;

          int noHits = inputCol->getNumberOfElements();
          int i_tmin; 
          std::vector<float>::iterator tmin;

          streamlog_out(DEBUG) << "Opened collection " << _inputColName << " contains " << noHits << " hits" << endl;

          for (int i = 0; i < noHits; i++){
            SimCalorimeterHit *aSimCalorimHit = dynamic_cast<SimCalorimeterHit*>(inputCol->getElementAt(i));

            int noSubHits = aSimCalorimHit->getNMCContributions();

            /* ---------------------------------- */

            vector<float> subhitsTime;
            vector<float> subhitsEnergy;
            subhitsTime.resize(noSubHits);
            subhitsEnergy.resize(noSubHits);

            // const int nImps = 4;
            vector<Imp> Imps(noSubHits);

            // Find out hit layer
            // vector<float>::iterator itr = std::find(_FixedPosZ_float.begin(), _FixedPosZ_float.end(), aSimCalorimHit->getPosition()[2]);
            // unsigned int hit_layer_index = std::distance(_FixedPosZ_float, itr);
            //int hit_layer_index = std::distance(_FixedPosZ_float.begin(), std::find(_FixedPosZ_float.begin(), _FixedPosZ_float.end(), aSimCalorimHit->getPosition()[2]));
            int hit_layer_index=9999;
	    vector<float> hitZtolayer;
	    float smallestdistance=9999.;
	    float hitZ=aSimCalorimHit->getPosition()[2];
	    int NLayers=0;
	    if(_FixedPosZ.size() != 0) { NLayers=_FixedPosZ.size(); }
	    else if(_NSlabs != 9999.) { NLayers=_NSlabs; }
	    else { streamlog_out(ERROR)<<"Missing number of slabs"<<std::endl; }
	    
	    for (int ilayer = 0; ilayer < NLayers; ilayer++){
	      hitZtolayer.push_back(abs(_FixedPosZ_float[ilayer]-hitZ));
	      //streamlog_out(DEBUG)<<"Layer comparing: "<<ilayer<<", Distance: "<<hitZtolayer.at(ilayer)<<endl;
	      if(ilayer==0) {
		smallestdistance=hitZtolayer.at(0);
		hit_layer_index=0;
	      }
	      else if(hitZtolayer.at(ilayer) == 0.) {
		smallestdistance=hitZtolayer.at(ilayer);
                hit_layer_index=ilayer;
		break;
	      }
	      else if(hitZtolayer.at(ilayer) < smallestdistance) {
		smallestdistance=hitZtolayer.at(ilayer);
		hit_layer_index=ilayer;
	      }
	    }
	    streamlog_out(DEBUG)<<"Closest slab: "<<hit_layer_index<<". Distance: "<<hitZtolayer.at(hit_layer_index)<<endl;
	    if(hitZtolayer.at(hit_layer_index) > _deltaZ) streamlog_out(WARNING)<<"Hit distance to layer over desired tolerance of "<<_deltaZ<<"mm"<<endl;
	      
	    float this_MIPvalue = _MIP2GeVFactors_float[hit_layer_index]; 
            if(evtNumber%1000 == 0) streamlog_out(DEBUG)<<"hit layer index: "<<hit_layer_index<<endl;
	    if(evtNumber%1000 == 0) streamlog_out(DEBUG)<<"Hit position: "<<aSimCalorimHit->getPosition()[2]<<endl;
            
	    for(int j = 0; j < noSubHits; j++){//fill hit time and deposited energy per sub-hit

              // const float *_hitStep = aSimCalorimHit->getStepPosition(j);

              // Float_t _hitOntileX =(Float_t)_hitStep[0]-(Float_t) std::floor(_hitStep[0]/_tileEdgeX)*_tileEdgeX;
              // Float_t _hitOntileY =(Float_t)_hitStep[1]-(Float_t) std::floor(_hitStep[1]/_tileEdgeY)*_tileEdgeY;

              float _edepstep = (float)aSimCalorimHit->getEnergyCont(j);
              // TODO: Check dead space!
              // if(_hitOntileX < halfdS || _hitOntileX > (_tileEdgeX-halfdS) || _hitOntileY < halfdS || _hitOntileY > (_tileEdgeY-halfdS)) _edepstep = 0.0;

              Imps[j] = {
                (float)aSimCalorimHit->getTimeCont(j),
                _edepstep/this_MIPvalue
	      };
	      //if(evtNumber%1000 == 0) cout<<"subhitenergy: "<<_edepstep<<", thisMIPvalue: "<<this_MIPvalue<<endl;
	      if(evtNumber%1000 == 0) streamlog_out(DEBUG)<< "\tsubhit " << j << ", energy/mip = " << (double)aSimCalorimHit->getEnergyCont(j) << "/" << this_MIPvalue << endl;
              subhitsTime[j] = (float)aSimCalorimHit->getTimeCont(j);
              subhitsEnergy[j] = _edepstep;

            }//end loop on subhits

            // Test ievent
            // if (_nEvt % 10 == 0) cout << "Passed loop subhit for evt " << _nEvt << endl;
            if ((_nEvt == 0) && (i == 0)) _plotHit  = true;
            else _plotHit = false;

            //PulseShapeCRRC(Imps, _plotHit);
            PulseShapeCRRC(Imps);
            _nTot++;
            if (_SSRead == 0.) {
              _nLeft++;
              continue;
            }
            // Plot shaper
            CalorimeterHitImpl * aCalorimHit =  new CalorimeterHitImpl();

            const float *hitpos = aSimCalorimHit->getPosition();

            aCalorimHit->setTime(_FSTime);
            aCalorimHit->setEnergy(_SSRead);
            aCalorimHit->setPosition(hitpos);
            aCalorimHit->setCellID0(aSimCalorimHit->getCellID0());

            calOutVec->addElement(aCalorimHit);


            if(evtNumber%1000 == 0) streamlog_out(DEBUG) << "Added hit to new collection " << _calorimOutCollections.at(icol) << endl;

            time_sha = _FSTime;
            energy_sha = _SSRead;
            posx = hitpos[0];
            posy = hitpos[1];
            posz = hitpos[2];


            // vector<int> lookup;
            // lookup.resize(noSubHits);

            // TMath::Sort(noSubHits,&subhitsTime[0],&lookup[0],kFALSE);//creating look-up table for time ordered subhits
            // time_sim = subhitsTime[lookup[0]];
            // energy_sim = subhitsEnergy[lookup[0]];
            //
            tmin = std::min_element(subhitsTime.begin(), subhitsTime.end()); 
            // i_tmin = TMath::LocMin(noSubHits, subhitsTime);
            i_tmin = std::distance(subhitsTime.begin(), tmin);
            time_sim = subhitsTime[i_tmin];
            //energy_sim = subhitsEnergy[i_tmin];
            // More reasonable to try with sum on all those arrays
            energy_sim = aSimCalorimHit->getEnergy();
            
            _treeout->Fill();
            


            // float energy(0.);//energy sum
            // int ThrIndex(0);//index of the first hit in the sliding window
            // float tH(0.);//time of passing threshold

            // float Epar(0.);
            // bool passThr(false);

            // /*checking for the time of signal passing threshold:
            // amplitude are added only until they are in a 50 ns window*/


            // for(int k=0; k<noSubHits; k++){//loop on all subhits
            //   ThrIndex = k;
            //   Epar = 0.;
            //   for(int isub=k; isub< noSubHits; isub++){//loop on the remaining hits

            //     if((subhitsTime[lookup[isub]]-subhitsTime[lookup[ThrIndex]]) <_tfast){//check if the following hit is within tfast
            //       Epar += subhitsEnergy[lookup[isub]];
            //     }else{break;}//endif //if time is > tfast exit loop on remaining hits

            //     if(Epar>energyThreshold){//check of energy sum vs threshold
            //       tH = subhitsTime[lookup[isub]];//ste time hit = last hit time
            //       passThr = true;
            //       break;//exit loop on remaining hits
            //     } //endif

            //   }//end loop on remaining hits
            //   if(passThr==true) break; //if threshold had passed exit loop on all subhits
            // }//end loop


            // /*checking if the subhits happen within the risetime of the slow shaper*/
            // if(passThr==true){
            //   for(int j=ThrIndex; j<noSubHits; j++){
            //     if(subhitsTime[lookup[j]] < (subhitsTime[lookup[ThrIndex]]+_tslow)){
            //       energy += subhitsEnergy[lookup[j]];
            //     }
            //   }

            //   CalorimeterHitImpl * aCalorimHit =  new CalorimeterHitImpl();

            //   const float *hitpos = aSimCalorimHit->getPosition();

            //   aCalorimHit->setTime(tH);
            //   aCalorimHit->setEnergy(energy);
            //   aCalorimHit->setPosition(hitpos);
            //   aCalorimHit->setCellID0(aSimCalorimHit->getCellID0());

            //   calOutVec->addElement(aCalorimHit);

            //   streamlog_out(DEBUG) << "Added hit to new collection " << _calorimOutCollections.at(icol) << endl;

            // }
          }//end loop over SimCalorimeterHits

          //=======================================================================

          //LCParameters &theParam = calOutVec->parameters();
          // theParam.setValue(LCIO::CellIDEncoding, _encoding);

	  //Insertion JP start
	  LCParameters &ogParam = inputCol->parameters();
          string InColEncoder=ogParam.getStringVal("CellIDEncoding");
          LCParameters &theParam = calOutVec->parameters();
          theParam.setValue("CellIDEncoding",InColEncoder);
	  //Insertion JP end

          evt->addCollection( calOutVec, _calorimOutCollections.at(icol) ) ;

          //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

          streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber()
          << "   in run:  " << evt->getRunNumber() << std::endl ;

          _nEvt ++ ;

        }//end if collection found
      }//end if find collection name
    }//end loop collections
  }

  /************************************************************************************/

  void ShapingProcessor::check( LCEvent * evt ) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
  }

  /************************************************************************************/
  void ShapingProcessor::end(){

    // Put here the linear fit of sim/shaped energy
    //EnergyLinearity();

    _treeout->Write();
    _rootout->Write("", TObject::kOverwrite);
    _rootout->Close();

    // cout << "Number of hits left out: " << to_string(_nLeft) << "\n";
    // cout << "Number of total hits: " << to_string(_nTot) << "\n";
    // cout << "Ratio: " << to_string((_nTot - _nLeft) / _nTot) << "\n";

    std::cout << "ShapingProcessor::end()  " << name()
    << " processed " << _nEvt << " events in " << _nRun << " runs "
    << std::endl ;
  }


  /************************************************************************************/
  void ShapingProcessor::printParameters(){
    std::cerr<<"============= Shaping Processor ================="<<std::endl;
    std::cerr<<" Simulating SiWECAL readout shaping and energy measurement (based on Ahc2)"<<std::endl;
    // std::cerr<<" Tile EdgeX: "<<_tileEdgeX <<" mm"<< " Tile EdgeY: "<<_tileEdgeY <<" mm" << std::endl;
    // std::cerr<<" Dead space bet. tiles: "<<_deadSpace<<" mm"<<std::endl;
    // std::cerr<<" disctiminator: "<<_MIPThr<<" MIP"<<std::endl;
    // std::cerr<<" fast shaper time: "<<_tfast<<" ns"<<std::endl;
    // std::cerr<<" slow shaper time: "<<_tslow<<" ns"<<std::endl;
    std::cerr<<"======================================================="<<std::endl;
    return;

  }

  /************************************************************************************/
  double ShapingProcessor::shaper(double x, double tau, int n, vector<Imp> Imps){
  
    tau = tau / n;
    double s = 0;
    for (auto& hit: Imps) { 
      double T = (x - hit.t) / tau;
      if (T > 0 && n > 0){
        s += 4 * hit.A / TMath::Factorial(n) * TMath::Power(T, n) * TMath::Exp(-T);
      }
    }
    return s;
  }

  /************************************************************************************/
  vector<struct ShapingProcessor::Imp> ShapingProcessor::histInput(vector<Imp> Imps){
    auto hInput = new TH1D("hImp", "Input hit histo", _nbinsI, 0, _nbinsI*_bwI);
    for (auto& hit: Imps) hInput->Fill(hit.t, hit.A);
    vector<Imp> hImps;
    for (int i = 0; i < _nbinsI; i++) {
      auto binContent = hInput->GetBinContent(i);
      if (binContent > 0) hImps.push_back({hInput->GetBinCenter(i), binContent});
    }
    
    delete hInput;
    return hImps;

  }

  
  /************************************************************************************/
  //void ShapingProcessor::PulseShapeCRRC(vector<Imp> Imps, bool _plotHit){
  void ShapingProcessor::PulseShapeCRRC(vector<Imp> Imps){
    
    // for (auto& hit: Imps) cout << "Hit hist:\t A = " << hit.A << "\tt = " << hit.t << "\n";
    
    if (_useHistInput) Imps = histInput(Imps);
    
    // for (auto& hit: Imps) cout << "Hit hist:\t A = " << hit.A << "\tt = " << hit.t << "\n";
    // if (Imps.size() >= 2) cout << "Hit hist of len " << Imps.size() << "\n";

    auto h_FS = new TH1D("h_FS", "FS Histo", _nbinsF, 0, _nbinsF*_bwF);
    auto h_SS = new TH1D("h_SS", "SS Histo", _nbinsS, 0, _nbinsS*_bwS);
    double shapedBin;

    for (int ibin = 0; ibin < _nbinsF; ibin++){
      shapedBin = shaper(h_FS->GetBinCenter(ibin), _tauF, _nF, Imps);
      if (_filterNoise) shapedBin = gRandom->Gaus(shapedBin, _FSNoise);
      h_FS->SetBinContent(ibin, shapedBin);

    }
    
    for (int ibin = 0; ibin < _nbinsS; ibin++){
      shapedBin = shaper(h_SS->GetBinCenter(ibin), _tauS, _nS, Imps);
      if (_filterNoise) shapedBin = gRandom->Gaus(shapedBin, _SSNoise);
      h_SS->SetBinContent(ibin, shapedBin);
    }

    if (h_FS->FindFirstBinAbove(_MIPThreshold) >= 0) {
      maxFS = h_FS->GetMaximum();
      maxSS = h_SS->GetMaximum();
      t_maxFS = h_FS->GetBinCenter(h_FS->GetMaximumBin());
      t_maxSS = h_SS->GetBinCenter(h_SS->GetMaximumBin());
    }

    else {
      maxFS = 0.; 
      maxSS = 0.;
      t_maxFS = 0.; 
      t_maxSS = 0.;
    }
    
    //for( int iThr = 0; iThr<_nThr; iThr++){
      //int BinThr = h_FS->FindFirstBinAbove(_Thrs[iThr]);
    int BinThr = h_FS->FindFirstBinAbove(_MIPThreshold);
    if (BinThr >= 0) {
      _FSTime = h_FS->GetBinCenter(BinThr);
      _SSRead = h_SS->GetBinContent(h_SS->FindBin(_FSTime + _delay));

    }
    else {
      _FSTime = 0.;
      _SSRead = 0.;
    }

    // }
    if (_plotHit and not(_plotted)) {
      TCanvas *c1 = new TCanvas("c1","Shaper Draw test", 600, 600);
      //if (_useHistInput) {
        auto hInput = new TH1D("hImp_plot", "Input hit histo", _nbinsI, 0, _nbinsI*_bwI);
        for (auto& hit: Imps) hInput->Fill(hit.t, hit.A);
        hInput->SetLineColor(kBlack);
        hInput->Draw("HIST");
      //}
      string text;
      text = "Input max: " + to_string(hInput->GetMaximum());

      h_FS->Draw("SAME");
      h_SS->Draw("SAME");
      hInput->GetXaxis()->SetRangeUser(0, 1.1*t_maxSS);
      hInput->GetYaxis()->SetRangeUser(0, 1.3*hInput->GetMaximum());
      c1->Update();
      TLine *line_thr = new TLine(0, _MIPThreshold, 1.1*t_maxSS, _MIPThreshold);
      line_thr->SetLineStyle(9);
      TLine *line_maxSS = new TLine(t_maxSS, 0, t_maxSS, maxSS);
      line_maxSS->SetLineStyle(9);
      TLine *line_read = new TLine(_FSTime, _SSRead, _FSTime+_delay, _SSRead);
      line_read->SetLineWidth(4);
      line_read->SetLineColor(kRed);
      TLine *line_toread = new TLine(_FSTime, _MIPThreshold, _FSTime, maxSS);
      line_toread->SetLineStyle(9);
      line_toread->SetLineColor(kRed);

      TPaveText *pt = new TPaveText(130, .3, 200, 1);
      pt->SetTextSize(0.03);
      pt->AddText(text.c_str());
      text = "Fast S. time: " + to_string(_FSTime);
      pt->AddText(text.c_str());
      text = "Slow S. energy: " + to_string(_SSRead);
      pt->AddText(text.c_str());
      //text = "Slow S energy: " + to_string(_SSRead) + "\n";
      //pt->AddText(text.c_str());
      pt->Draw();


      line_thr->Draw("SAME");
      line_maxSS->Draw("SAME");
      line_read->Draw("SAME");
      line_toread->Draw("SAME");
      c1->SaveAs("hit_shaper.root");
      _plotted = true;
    }

  delete h_FS;
  delete h_SS;

  }
  // void ShapingProcessor::EnergyLinearity() {
  // 
  // TTreeReader reader(_treeout);
  // TTreeReaderValue<double> energy_sim(reader, "energy_sim");
  // TTreeReaderValue<double> energy_sha(reader, "energy_sha");

  // int i = 0;
  // auto N = reader.GetEntries();
  // double x[N], y[N];
  // while (reader.Next() and i < N-1) {
  // //while (reader.Next()) {
  //   //std::cout << *energy_sim << " ";
  //   x[i]=*energy_sim;
  //   y[i]=*energy_sha;
  //   i++;
  //   } 
  // TF1 *f1 = new TF1("f1","[0]+[1]*x");
  // TGraph *gr = new TGraph(N, x, y);
  // gr->Fit("f1");
  // std::cout << "Params " << f1->GetParameter(1);
  // 
  // }

}
