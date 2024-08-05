#include "ShapedConversionProcessor.hh"

#include <iostream>
#include <cstdlib>
#include <unistd.h>


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
#include <TF1Convolution.h>
//#include <TString.h>
#include <TFile.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TDirectory.h>
#include <TSpectrum.h>
#include <TRandom.h>
#include "TMath.h"
#include "TROOT.h"


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

  ShapedConversionProcessor aShapedConversionProcessor;

  ShapedConversionProcessor::ShapedConversionProcessor() : Processor("ShapedConversionProcessor")
  {

    vector<string> calorimInpCollectionsExample;
    registerProcessorParameter("Input_Collections",
                               "vector of collections for input",
                               _calorimInpCollections,
                               calorimInpCollectionsExample);
    
    vector<string> siThicknessesExample;
    registerProcessorParameter("SiThicknesses",
                               "vector of Silicon thicknesses per layer",
                               _siThicknesses,
                               siThicknessesExample);


    string eConfNameExample;
    registerProcessorParameter("Energy_Conf_Name",
                               "Name to identify auxiliary file",
                               _eConfName,
                               eConfNameExample);
    
    registerProcessorParameter("MIPFitMode",
                               "Fit mode of MIP in cells (1, 2, or 3)",
                               _MIPFitMode,
                               int(1));
    
    
  }

  /************************************************************************************/

  void ShapedConversionProcessor::init()
  {
    streamlog_out(DEBUG) << "init called" << std::endl ;

    //K Map
    //_mapKs.clear();
    //FillMapK();
    printParameters();

    _nRun = 0 ;
    _nEvt = 0 ;
    _maxSubHits = 0;
    _maxHits = 0;
    _writeHisto_mu_thr = 0.000098; // May be a steering file? Or maybe not if params are stored in ttrees
    _writeHisto = false;
    // Wrong tiling
    // for (float val = -88.0; val <= 88.0; val += 5.5) _FixedPosXY.push_back(val);
    for (float val = -85.25; val <= 85.25; val += 5.5) _FixedPosXY.push_back(val);
    // Same for the three confs
    // _FixedPosZ = {12.327500, 27.327499, 42.327499, 57.327499,
    //               72.327499, 87.327499, 147.327499};
    _FixedPosZ = {6.225,  21.225,  36.15,  51.15,  66.06,  81.06,  96.06,\
                  111.15, 126.15, 141.15, 156.15, 171.06, 186.06, 201.06,\
                  216.06};
    
    // Define unique thickness values (2017: single value, 2021: two, 2022: three)
    // vector<float> float_siThicknesses(_siThicknesses.begin(), _siThicknesses.end());
    vector<float> float_siThicknesses(_FixedPosZ.size());
    for (int iz = 0; iz < _FixedPosZ.size(); iz++)
    {
      float_siThicknesses[iz] = stof(_siThicknesses[iz]);
    }
    auto _unique_siThicknesses = float_siThicknesses;
    auto it = unique(float_siThicknesses.begin(), float_siThicknesses.end());
    _unique_siThicknesses.resize(distance(float_siThicknesses.begin(), it));
    

    // Prototype for TTree in params
    // TTree *_fitParams = new TTree("FitParams","Energy cell fit parameters");
    // _gauss_c =      _fitParams->Branch("_gauss_c",&_gauss_c,"_gauss_c/D"); 
    // _gauss_mu =     _fitParams->Branch("_gauss_mu",&_gauss_mu,"_gauss_mu/D"); 
    // _gauss_sigma =  _fitParams->Branch("_gauss_sigma",&_gauss_sigma,"_gauss_sigma/D"); 
    // _landau_c =     _fitParams->Branch("_landau_c",&_landau_c,"_landau_c/D"); 
    // _landau_mpv =   _fitParams->Branch("_landau_mpv",&_landau_mpv,"_landau_mpv/D"); 
    // _landau_sigma = _fitParams->Branch("_landau_sigma",&_landau_sigma,"_landau_sigma/D"); 
    // _langaus_c =    _fitParams->Branch("_langaus_c",&_langaus_c,"_langaus_c/D"); 
    // _chi2ndf =      _fitParams->Branch("_chi2ndf",&_chi2ndf,"_chi2ndf/D"); 

    
    // Fit histos, gauss, landau
    _gauss_c = TH1F("gauss_c", "Gauss constant distribution", 100, 0., 6000.);
    _gauss_mu = TH1F("gauss_mu", "Gauss mu distribution", 400, 0.05, 10.);
    _gauss_sigma = TH1F("gauss_sigma", "Gauss sigma distribution", 100, 0., 3.);
    _landau_c = TH1F("landau_c", "Landau constant distribution", 100, 0., 6000.);
    // For W confs
    //_landau_mpv = TH1F("landau_mpv", "Landau MPV distribution", 100, 0.00005, 0.00025);
    // for conf0
    _landau_mpv = TH1F("landau_mpv", "Landau MPV distribution", 400, 0.05, 10.);
    _landau_mpv_cor = TH1F("landau_mpv_cor", "Landau Corrected MPV distribution", 400, 0., 10.);
    _landau_sigma = TH1F("landau_sigma", "Landau sigma distribution", 100, 0., 3.);
    _chi2ndf = TH1F("chi2ndf", "chi2ndf distribution", 200, 0., 20);
    // Langaus
    _langaus_c = TH1F("langaus_c", "Landau (x) Gaussian constant distribution", 100, 0., 6000.);


    for (unsigned int ix = 0; ix < _FixedPosXY.size(); ix++)
    {
      for (unsigned int iy = 0; iy < _FixedPosXY.size(); iy++)
      {
        for (unsigned int iz = 0; iz < _FixedPosZ.size(); iz++)
        {
          std::string ename = "Energy_GeV__" + std::to_string(ix) + "_"\
                                             + std::to_string(iy) + "_"\
                                             + std::to_string(iz);
          energyHistos[ix][iy][iz] = new TH1F(ename.c_str(), ename.c_str(),
                                              //100, 0., 0.003);
                                              //200, 0., 0.001); // default for long time
                                              100, 0., 10.);
          std::string tname = "Time_ns__" + std::to_string(ix) + "_"\
                                          + std::to_string(iy) + "_"\
                                          + std::to_string(iz);
          _timeHistos[ix][iy][iz] = new TH1F(tname.c_str(), tname.c_str(),
                                              1600, 1., 5);
        }
      }
    }
    for (unsigned int iz = 0; iz < _FixedPosZ.size(); iz++)
    {
      std::string lname = "Energy_GeV__layer_" + std::to_string(iz);
      layerEnergyHistos[iz] = new TH1F(lname.c_str(), lname.c_str(),
                                          100, 0., 10.);

      std::string tname_layer = "Time_ns_layer_" + std::to_string(iz);
      _timeHistos_layer[iz] = new TH1F(tname_layer.c_str(), tname_layer.c_str(),
                                           1600, 1., 5);
    }
    //_rootout = new TFile("test_energies.root", "RECREATE");
    _rootout = new TFile(_eConfName.c_str(), "RECREATE");
    _high_mu_dir = _rootout->mkdir("high_mu");
    
    // int randint;
    // _randIndices.push_back(0);
    // gRandom = new TRandom(4357);
    // while (_randIndices.size() <= 100)  
    // {
    //   randint = gRandom->Integer(1000);
    //   if (std::find(_randIndices.begin(), _randIndices.end(), randint) == _randIndices.end())
    //     _randIndices.push_back(randint);
    // }

  }


  /************************************************************************************/
  void ShapedConversionProcessor::processRunHeader( LCRunHeader* run)
  {
    _nRun++ ;
  }

  void ShapedConversionProcessor::processEvent( LCEvent * evt )
  {

    // Fabricio
    //std::vector<int> CellIDs;


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

          // LCCollectionVec* calOutVec = new LCCollectionVec( LCIO::SIMCALORIMETERHIT) ;

          // LCFlagImpl hitFlag(calOutVec->getFlag());
          // hitFlag.setBit(LCIO::RCHBIT_TIME);
          // hitFlag.setBit(LCIO::CHBIT_STEP);
          // hitFlag.setBit(LCIO::CHBIT_LONG);
          // calOutVec->setFlag(hitFlag.getFlag());

          int noHits = inputCalorimCollection->getNumberOfElements();
          //if (_maxHits < noHits) _maxHits = noHits;
          //_encoding = inputCalorimCollection->getParameters().getStringVal("CellIDEncoding");

          //CellIDDecoder<SimCalorimeterHit> decoder(inputCalorimCollection);
          //CellIDEncoder<SimCalorimeterHitImpl> encoder(_encoding.c_str(), calOutVec);

          // Likely useless for our purposes
          // bool hasKminus1 = false;
          // if (_encoding.find("K-1") != std::string::npos)
          // {
          //   hasKminus1 = true;
          // }

          for (int i = 0; i < noHits; i++)
          {
            SimCalorimeterHit *aSimCalorimHit = dynamic_cast<SimCalorimeterHit*>(inputCalorimCollection->getElementAt(i));
            //New SimHit converted to MIP
            // Fabricio: line below useless?
            SimCalorimeterHitImpl *aSimCalorimHit_MIP =  new SimCalorimeterHitImpl();

            int noSubHits = aSimCalorimHit->getNMCContributions();
            if (_maxSubHits < noSubHits) _maxSubHits = noSubHits;
            for(int j = 0; j < noSubHits; j++)
            {
              EVENT::MCParticle *mcp_subhit = (EVENT::MCParticle*)aSimCalorimHit->getParticleCont(j);
              float *step_subhit = (float*)aSimCalorimHit->getStepPosition(j);
              int pdg_subhit = (int)aSimCalorimHit->getPDGCont(j);
              float time_subhit = (float)aSimCalorimHit->getTimeCont(j);
              float edep_subhit = (float)aSimCalorimHit->getEnergyCont(j);

              //Adding contribution to new hit
              aSimCalorimHit_MIP->addMCParticleContribution(mcp_subhit, edep_subhit/_MIPvalue, time_subhit, pdg_subhit, step_subhit);
            }
            // Deleted commented encoder/decoder I/J/K stuff


            // Fabricio
            int this_cellid = aSimCalorimHit->getCellID0();
            float this_cellPosX = aSimCalorimHit->getPosition()[0];
            float this_cellPosY = aSimCalorimHit->getPosition()[1];
            float this_cellPosZ = aSimCalorimHit->getPosition()[2];
            if (std::find(_CellIDs.begin(), _CellIDs.end(), this_cellid) == _CellIDs.end())
              _CellIDs.push_back(this_cellid);

            if (std::find(_CellPosX.begin(), _CellPosX.end(), this_cellPosX) == _CellPosX.end())
              _CellPosX.push_back(this_cellPosX);
            if (std::find(_CellPosY.begin(), _CellPosY.end(), this_cellPosY) == _CellPosY.end())
              _CellPosY.push_back(this_cellPosY);
            if (std::find(_CellPosZ.begin(), _CellPosZ.end(), this_cellPosZ) == _CellPosZ.end())
              _CellPosZ.push_back(this_cellPosZ);

            int ix = 0;
            int iy = 0;
            int iz = 0;

            std::vector<float>::iterator itrx;
            itrx = std::find(_CellPosX.begin(), _CellPosX.end(), this_cellPosX);
            if (itrx != _CellPosX.cend()) ix = std::distance(_CellPosX.begin(), itrx);
            else std::cout << "Cell X not found";
           
            std::vector<float>::iterator itry;
            itry = std::find(_CellPosY.begin(), _CellPosY.end(), this_cellPosY);
            if (itry != _CellPosY.cend()) iy = std::distance(_CellPosY.begin(), itry);
            else std::cout << "Cell Y not found";
           
            std::vector<float>::iterator itrz;
            itrz = std::find(_CellPosZ.begin(), _CellPosZ.end(), this_cellPosZ);
            if (itrz != _CellPosZ.cend()) iz = std::distance(_CellPosZ.begin(), itrz);
            else std::cout << "Cell Z not found";

            energyHistos[ix][iy][iz]->Fill(aSimCalorimHit->getEnergy());

            for (int it = 0; it < noSubHits; it++)
            {
              _timeHistos[ix][iy][iz]->Fill(aSimCalorimHit->getTimeCont(it),
                                            aSimCalorimHit->getEnergyCont(it));
              _timeHistos_layer[iz]->Fill(aSimCalorimHit->getTimeCont(it),
                                          aSimCalorimHit->getEnergyCont(it));
            }


            //Setting new Hit
            
            // aSimCalorimHit_MIP->setCellID0(aSimCalorimHit->getCellID0());
            // //aSimCalorimHit_MIP->setEnergy(aSimCalorimHit->getEnergy()/_MIPvalue);
            // aSimCalorimHit_MIP->setEnergy(aSimCalorimHit->getEnergy());
            // aSimCalorimHit_MIP->setPosition(aSimCalorimHit->getPosition());

            // calOutVec->addElement(aSimCalorimHit_MIP);

          }//end loop over SimCalorimeterHits

          // LCParameters &theParam = calOutVec->parameters();
          // theParam.setValue(LCIO::CellIDEncoding, _encoding);

          //if(calOutVec->getNumberOfElements() > 0)
          //evt->addCollection( calOutVec, _calorimOutCollections.at(icol) ) ;

          //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

          streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber()
          << "   in run:  " << evt->getRunNumber() << std::endl ;

          _nEvt ++ ;

          if (_nEvt % 1000 == 0) _rootout->Write("", TObject::kOverwrite);
        }//end if col
      }//end find col names
    }//end for loop
    std::sort(_CellPosX.begin(), _CellPosX.end());
    std::sort(_CellPosY.begin(), _CellPosY.end());
    std::sort(_CellPosZ.begin(), _CellPosZ.end());
  }

  /************************************************************************************/

  void ShapedConversionProcessor::check( LCEvent * evt ) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
  }

  /************************************************************************************/
  void ShapedConversionProcessor::end(){

    // Setting fit range and start values
    // double fr[2];
    // double sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    // double chisqr;
    // int    ndf;
    
    // int maxBin;
    // int shift = 2;
    int writeCounter = 0;
    // double minFit, maxFit;

    for (unsigned int iz = 0; iz < _FixedPosZ.size(); iz++)
    {
       for (unsigned int ix = 0; ix < _FixedPosXY.size(); ix++)
       {
         for (unsigned int iy = 0; iy < _FixedPosXY.size(); iy++)
         {
           layerEnergyHistos[iz]->Add(energyHistos[ix][iy][iz]);
           delete energyHistos[ix][iy][iz];
         }
       }
       // for (unsigned int ix = 0; ix < _FixedPosXY.size(); ix++)
       // {
       //   for (unsigned int iy = 0; iy < _FixedPosXY.size(); iy++)
       //   {
       //     if (energyHistos[ix][iy][iz]->GetEntries() > 500)
       //     {
            
             // if (_MIPFitMode == 1) {fitGaus(energyHistos[ix][iy][iz]);}
             // else if (_MIPFitMode == 2) {fitLandau(energyHistos[ix][iy][iz]);}
             // else if (_MIPFitMode == 3) {fitLanGaus(energyHistos[ix][iy][iz]);}
             // else
             // {
             //   cout << "Invalid _MIPFitMode = " << std::to_string(_MIPFitMode) << "\n";
             //   break;
             // }
             
             if (_MIPFitMode == 1) {fitGaus(layerEnergyHistos[iz]);}
             else if (_MIPFitMode == 2) {fitLandau(layerEnergyHistos[iz]);}
             else if (_MIPFitMode == 3) {fitLanGaus(layerEnergyHistos[iz]);}
             else
             {
               cout << "Invalid _MIPFitMode = " << std::to_string(_MIPFitMode) << "\n";
               break;
               }

             layerEnergyHistos[iz]->Write();
             // if (writeCounter < 20 && !_writeHisto)
             // {
             //   layerEnergyHistos[iz]->Write();
             //   writeCounter++;
             //  
             // }
             // if (_writeHisto)
             // {
             //   // TDirectory *_high_mu_dir = top->mkdir("high_mu");
             //   _high_mu_dir->cd();
             //   energyHistos[ix][iy][iz]->Write();
             //   _rootout->cd();
             // }

             //_timeHistos[ix][iy][iz]->Write(); // atm per-cell timing histos not needed...

           // }
           delete layerEnergyHistos[iz];
       //   }
       // }
     //_timeHistos_layer[iz]->Write();
     }       
    for (unsigned int iz = 0; iz < _FixedPosZ.size(); iz++)
    {
     _timeHistos_layer[iz]->Write();
    }

    // Fit box:
    //auto style = gROOT->GetStyle(style_name);

    _gauss_c.Write();
    _gauss_mu.Write();
    _gauss_sigma.Write();
    _landau_c.Write();
    _landau_mpv.Write();
    _landau_mpv_cor.Write();
    _landau_sigma.Write();
    _langaus_c.Write();
    _chi2ndf.Write();

    _rootout->Write("", TObject::kOverwrite);
    _rootout->Close();
    //    
    //    std::cout << "Max hits in an event, max subhits in a hit " << _maxHits << ", " << _maxSubHits;
    
    
    std::cout << "ShapedConversionProcessorProcessor::end()  " << name()
      << " processed " << _nEvt << " events in " << _nRun << " runs "
      << "FLAG"
      << std::endl ;

  }


  /************************************************************************************/
  void ShapedConversionProcessor::printParameters(){
    std::cerr << "============= SiWECALGeV2MIP Processor =================" <<std::endl;
    //std::cerr << "Converting Simulation Hits from GeV to MIP" <<std::endl;
    std::cerr << "Converting Simulation Hits from GeV to MIP (SiWECAL test)" <<std::endl;
    std::cerr << "Input Collection name : "; for(unsigned int i = 0; i < _calorimInpCollections.size(); i++) std::cerr << _calorimInpCollections.at(i) << " ";
    std::cerr << std::endl;
    //std::cerr << "Output Collection name : "; for(unsigned int i = 0; i < _calorimOutCollections.size(); i++) std::cerr << _calorimOutCollections.at(i) << " ";
    std::cerr << std::endl;
    std::cerr << "=======================================================" <<std::endl;
    return;

  }

  void ShapedConversionProcessor::fitGaus(TH1F* h)
  {
    // int shift = 2; // Have this as steering?
    _writeHisto = false;
    int maxBin;
    double halfMax, minFit = 0., maxFit = 0.;
    maxBin = h->GetMaximumBin();
    halfMax = h->GetBinContent(maxBin)/2.;
    for (int ibin = 0; ibin < maxBin; ibin++)
    {
      if (h->GetBinContent(ibin) > halfMax)
      {
        minFit = h->GetBinCenter(ibin);
        break;
      }
    }

    for (int ibin = maxBin; ibin < h->GetNbinsX(); ibin++)
    {
      if (h->GetBinContent(ibin) < halfMax)
      {
        maxFit = h->GetBinCenter(ibin);
        break;
      }
    }

    // Syntax below (4 params, incl/ bin range) works in newer versions of ROOT
    // minFit = h->GetBinCenter(h->FindFirstBinAbove(halfMax, 1, 1, maxBin));
    // maxFit = h->GetBinCenter(h->FindLastBinAbove(halfMax, 1, maxBin, -1));
    
    // Older version
    // minFit = h->GetBinCenter(maxBin-shift);
    // maxFit = h->GetBinCenter(maxBin+shift);
    
    TF1 *f = new TF1("f", "gaus", minFit, maxFit);
    f->SetParameters(h->GetMaximum(), h->GetBinCenter(maxBin), (maxFit-minFit)/2);
    h->Fit("f", "RQ");
    gStyle->SetOptFit(1111);
    f->Draw("SAME");
    
    _gauss_c.Fill(f->GetParameter(0));
    _gauss_mu.Fill(f->GetParameter(1));
    if (f->GetParameter(1) > _writeHisto_mu_thr) {_writeHisto = true;}
    _gauss_sigma.Fill(f->GetParameter(2));
    _chi2ndf.Fill(f->GetChisquare()/f->GetNDF());
    cout << f->GetParameter(1) << " +- " << f->GetParameter(2) << endl;
  }

  void ShapedConversionProcessor::fitLandau(TH1F* h)
  {
    _writeHisto = false;
    double par[3], par_lo[3], par_hi[3];
    // Parameters
    par[0] = h->Integral()/10.;
    par[1] = h->GetMean();
    par[2] = h->GetRMS()/10.;
    // Limits
    par_lo[0] = h->Integral("width")/10.;
    par_hi[0] = h->GetEntries()*2;
    par_lo[1] = 0.00005;
    par_hi[1] = h->GetMean()*2;
    par_lo[2] = h->GetRMS()/100.;
    par_hi[2] = h->GetRMS();
    // Fit range
    double minFit, maxFit;
    minFit = TMath::Max(h->GetMean() - h->GetRMS(), 0.);
    maxFit = h->GetMean() + 2. * h->GetRMS();
    
    TF1 *f = new TF1("f", "[0]*TMath::Landau(x,[1],[2],1)", minFit, maxFit);
    // Set
    for(int i = 0; i < 3; i++)
    {
      f->SetParameter(i, par[i]);
      f->SetParLimits(i, par_lo[i], par_hi[i]);
    }
    // Fit
    h->Fit("f","RBQM");
    
    gStyle->SetOptFit(1111);
    f->Draw("SAME");
    
    _landau_c.Fill(f->GetParameter(0));
    _landau_mpv.Fill(f->GetParameter(1));
    double mpv_corrected = f->GetParameter(1) + 0.22278298 * f->GetParameter(2);
    _landau_mpv_cor.Fill(mpv_corrected);
    if (mpv_corrected > _writeHisto_mu_thr) {_writeHisto = true;}
    _landau_sigma.Fill(f->GetParameter(2));
    _chi2ndf.Fill(f->GetChisquare()/f->GetNDF());
  }
  
  // Landau (x) Gaus by Adrián - Likely from
  // https://accserv.lepp.cornell.edu/svn/packages/root/tutorials/fit/langaus.C

  // TODO: Should be renamed to comply w/ code standards
  double ShapedConversionProcessor::langaufun(double *x, double *par)
  {
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation), 
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.
    
    // Numeric constants
    double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    double mpshift  = -0.22278298;       // Landau maximum location
    
    // Control constants
    double np = 100.0;      // number of convolution steps
    double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
    
    // Variables
    double xx, mpc, fland, xlow, xupp, step, i;
    double sum = 0.0;
    
    
    // MP shift correction
    mpc = par[1] - mpshift * par[0]; 
    
    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    
    step = (xupp-xlow) / np;
    
    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
      
      xx = xupp - (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }
    
    return (par[2] * step * sum * invsq2pi / par[3]);
  }



  TF1* ShapedConversionProcessor::langaufit(TH1F *his, double *fitrange, double *startvalues, double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors, double *ChiSqr, int *NDF)
  
  {
     // Once again, here are the Landau * Gaussian parameters:
     //   par[0]=Width (scale) parameter of Landau density
     //   par[1]=Most Probable (MP, location) parameter of Landau density
     //   par[2]=Total area (integral -inf to inf, normalization constant)
     //   par[3]=Width (sigma) of convoluted Gaussian function
     //
     // Variables for langaufit call:
     //   his             histogram to fit
     //   fitrange[2]     lo and hi boundaries of fit range
     //   startvalues[4]  reasonable start values for the fit
     //   parlimitslo[4]  lower parameter limits
     //   parlimitshi[4]  upper parameter limits
     //   fitparams[4]    returns the final fit parameters
     //   fiterrors[4]    returns the final fit errors
     //   ChiSqr          returns the chi square
     //   NDF             returns ndf
  
     int i;
     Char_t FunName[100];
  
     sprintf(FunName,"Fitfcn_%s",his->GetName());
  
     TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
     if (ffitold) delete ffitold;
  
     TF1 *ffit = new TF1(FunName, langaufun, fitrange[0], fitrange[1], 4);
     ffit->SetParameters(startvalues);
     ffit->SetParNames("Width","MP","Area","GSigma");
     
     for (i=0; i<4; i++) {
        ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
     }
  
     his->Fit(FunName,"RBQM");   // fit within specified range, use ParLimits, quiet, improve fit results (TMINUIT)
  
     ffit->GetParameters(fitparams);    // obtain fit parameters
     for (i=0; i<4; i++) {
        fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
     }
     ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
     NDF[0] = ffit->GetNDF();           // obtain ndf
  
     return (ffit);              // return fit function
  
  }
  
  void ShapedConversionProcessor::fitLanGaus(TH1F* h)
  {

    _writeHisto = false;
    
    // Fromn https://github.com/SiWECAL-TestBeam/SiWECAL-TB-analysis/blob/TB201706/singleslab/singleSlabAnalysis.cc
    Double_t fr[2];
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];   
    
    // Width // MPV // Integral // GaussSigma
    pllo[0]=h->GetRMS()/20.;     pllo[1]=0.0;             pllo[2]=0.005; pllo[3]=0.;
    plhi[0]=h->GetRMS();         plhi[1]=h->GetMean()*2.; plhi[2]=3.;    plhi[3]=h->GetRMS();
    
	fr[0]=h->GetMean()-0.8*h->GetRMS();
    fr[1]=1.5*h->GetMean();//+0.1*mip_histo.at(ichip).at(ichn)->GetRMS();
    
	sv[0]=h->GetRMS()*0.5;
    sv[1]=h->GetMean()*0.6;
    sv[2]=h->Integral("width");
    sv[3]=h->GetRMS()/5.;


    //  // Setting fit range and start values
    //  double fr[2];
    //  double sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    //  //limits width
    //  pllo[0] = h->GetRMS()/100.; 
    //  plhi[0] = h->GetRMS();
    //  //limits MPV
    //  pllo[1] = 0.00005;
    //  plhi[1] = h->GetMean()*2.;
    //  //limits Integral
    //  pllo[2] = h->Integral("width")/10.; 
    //  plhi[2] = h->GetEntries()*2; 
    //  //limits gaussian width
    //  // pllo[3] = h->GetRMS()/10.;
    //  // plhi[3] = h->GetRMS()*10;
    //  pllo[3] = h->GetRMS()/20.;
    //  plhi[3] = h->GetRMS()*5;

    //  //range of fit
    //  fr[0] = TMath::Max(h->GetMean() - h->GetRMS(), 0.);
    //  fr[1] = h->GetMean() + 2. * h->GetRMS();

    //  //initial values
    //  sv[0] = h->GetRMS() / 10.;//initial value landau width
    //  sv[1] = h->GetMean(); //mpv
    //  sv[2] = h->Integral();
    //  sv[3] = h->GetRMS();
    double chisqr;
    Int_t    ndf;
    // double par[3], par_lo[3], par_hi[3];
    // // Parameters
    // par[0] = h->Integral()/10.;
    // par[1] = h->GetMean();
    // par[2] = h->GetRMS()/10.;
    // par[3] = h->GetMean();
    // // Limits
    // par_lo[0] = h->Integral("width")/10.;
    // par_hi[0] = h->GetEntries()*2;
    // par_lo[1] = 0.00005;
    // par_hi[1] = h->GetMean()*2;
    // par_lo[2] = h->GetRMS()/100.;
    // par_hi[2] = h->GetRMS();
    // par_lo[3] = 0.00005;
    // par_hi[3] = h->GetMean();
    // // Fit range
    // double minFit, maxFit;
    // minFit = TMath::Max(h->GetMean() - h->GetRMS(), 0.);
    // maxFit = h->GetMean() + 2. * h->GetRMS();
    // 
    // TF1 *landau = new TF1("landau", "[0]*TMath::Landau(x,[1],[2],1)", minFit, maxFit);
    // // TF1 *gaus = new TF1("gaus", "gaus", minFit, maxFit);
    // TF1 *gaus = new TF1("gaus", "gaus", minFit, maxFit);
    // TF1Convolution *f_conv = new TF1Convolution(landau, gaus, minFit, maxFit);
    // f_conv->SetNofPointsFFT(1000);
    // TF1 *f = new TF1("f", *f_conv, minFit, maxFit, f_conv->GetNpar());
    // // Set
    // for(int i = 0; i < 4; i++)
    // {
    //   f->SetParameter(i, par[i]);
    //   f->SetParLimits(i, par_lo[i], par_hi[i]);
    // }
    // // Fit
    // h->Fit("f","RBQM");
    // 
    // gStyle->SetOptFit(1111);
    // f->Draw("SAME");
    //

    TF1 *f = langaufit(h, fr, sv, pllo, plhi, fp, fpe, &chisqr, &ndf);
    gStyle->SetOptFit(1111);
    h->Draw();
    
    _langaus_c.Fill(f->GetParameter(2));
    _landau_mpv.Fill(f->GetParameter(1));
    double mpv_corrected = f->GetParameter(1) + 0.22278298 * f->GetParameter(0);
    _landau_mpv_cor.Fill(mpv_corrected);
    if (mpv_corrected > _writeHisto_mu_thr) {_writeHisto = true;}
    _landau_sigma.Fill(f->GetParameter(0));
    _gauss_sigma.Fill(f->GetParameter(3));
    // _chi2ndf.Fill(f->GetChisquare()/f->GetNDF());
    _chi2ndf.Fill(chisqr / ndf);

  }


}
