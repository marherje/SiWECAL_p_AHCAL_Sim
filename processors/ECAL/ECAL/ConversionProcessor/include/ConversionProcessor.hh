#ifndef ConversionProcessor_HH
#define ConversionProcessor_HH 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <map>

//ROOT
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TDirectory.h>
#include <TFile.h>

using namespace lcio ;
using namespace marlin ;
using std::vector;

namespace CALICE
{
  /**
  * @brief Class doing the conversion between GeV and MIP + Correction of K layer (not smart!!)
  *
  *
  * @author eldwan.brianne@desy.de
  * @version 1.0
  * @date November 2015
  */

  class ConversionProcessor : public Processor {

  private:
    //processor parameters:
    float _MIPvalue; /**<MIP conversion factor in GeV*/
    std::string _rootfilename; /**<rootfile name*/
    std::string _auxFileName; /**<econf name*/

  public:

    virtual Processor*  newProcessor() { return new ConversionProcessor ; }

    /**Constructor
    */
    ConversionProcessor() ;

    /**Destructor
    */
    ~ConversionProcessor() {};

    /**Initialise
    */
    virtual void init() ;

    /**Process header
    */
    virtual void processRunHeader( LCRunHeader* run ) ;

    /**Process event (the working horse)
    */
    virtual void processEvent( LCEvent * evt ) ;

    /**Check
    */
    virtual void check( LCEvent * evt ) ;

    /**End of processing
    */
    virtual void end() ;

    /**Print Parameters
    */
    virtual void printParameters();
    
    /* Fit functions */
    virtual float fitGaus(TH1F* h);
    virtual float fitLandau(TH1F* h);
    virtual float fitLanGaus(TH1F* h);
    //Langaus functions below should be renamed/standardized
    static double langaufun(double *x, double *par);
    virtual TF1* langaufit(TH1F *his, double *fitrange, double *startvalues, double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors, double *ChiSqr, int *NDF);
    /** Fill the K map
     */
    // virtual void FillMapK();


  protected:
    std::string _calorimInpCollection; /**<input collection name*/
    std::string _calorimOutCollection; /**<output collection name*/

    vector<std::string> _calorimInpCollections;/**<input collection name*/
    vector<std::string> _calorimOutCollections;/**<output collection name*/
    //vector<std::string> _siThicknesses; WHY? (JESUS)

    std::string _encoding;

    std::string _fileMapK; /** File containing the relation between K in Simulation and K in data*/
    // vector<std::string> _filepaths; /** File containing the relation between K in Simulation and K in data*/

    std::map<int,std::map<int,int> > _mapK; /**<map containing SimK and K to convert*/
    // std::vector< std::map<int, int> > _mapKs; /**<map containing SimK and K to convert*/

    int _nRun ; /**<run number*/
    int _nEvt ; /**<evt number*/
    int _maxHits, _maxSubHits; 
    int _MIPFitMode;
    double _writeHisto_mu_thr;
    bool _writeHisto;
    std::vector<int> _randIndices;
    std::vector<int> _CellIDs ;
    std::vector<float> _CellPosX ;
    std::vector<float> _CellPosY ;
    std::vector<float> _CellPosZ ;
    
    // Should probably be arrays
    std::vector<float> _FixedPosXY ;
    
    // New fixed layers input
    int _NSlabs;
    
    // For output
    std::vector<float> _FixedPosZ_float ;
    std::vector<float> _MIP2GeVFactors_float;
    
    // Forgot to prepstringend _, number of cells, layers should be set in steering
    TH1F* energyHistos[33][33][15];
    TH1F* layerEnergyHistos[15];
    TH1F* _timeHistos[33][33][15];
    TH1F* _timeHistos_layer[15];
    TFile* _rootout;
    TDirectory* _high_mu_dir;
    // Prototype for having those as ttrees
    // TTree* _fitParams;
    // double _gauss_c, _gauss_mu, _gauss_sigma;
    // double _landau_c, _landau_mu, _landau_sigma;
    // double _langaus_c;
    // double _chi2ndf;
    // Histogram implementation
    TH1F _gauss_c, _gauss_mu, _gauss_sigma;
    TH1F _landau_c, _landau_mpv, _landau_mpv_cor, _landau_sigma;
    TH1F _langaus_c;
    TH1F _chi2ndf;
  } ;
}
#endif
