#ifndef ShapingProcessor_hh
#define ShapingProcessor_hh 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

#include <TH1.h>
#include <TFile.h>
#include <TTree.h>

using namespace lcio ;
using namespace marlin ;
using namespace std;

namespace CALICE
{

  /**
  * @brief Class for doing the AHCAL ROC Simulating behavior
  *
  * @author eldwan.brianne@desy.de
  * @version 2.0
  * @date November 2015
  */

  class ShapingProcessor : public Processor
  {
  private:
    // //processor parameters:
    // float _tileEdgeX;/**<Tile size in X*/
    // float _tileEdgeY;/**<Tile size in Y*/
    // float _deadSpace;/**<dead space*/
    // float _MIPvalue;/**<MIP Value*/
    // float _MIPThr;/**<MIP Threshold*/
    // float _tfast;/**<fast shaper integration time*/
    // float _tslow;/**<slow shaper integration time*/
    // int _NLayer;/**<Number of layers*/
    // string _LayerPattern;/**<layer pattern*/
    // vector<int> _LayerPatternVector;/**<layer pattern vector*/
    
    // Fabricio
    //double _bwI, _tauF, _tauS, _delay, _bwF, _bwS;
    //int _nbinsI, _nbinsS, _nbinsF, _nF, _nS;
    //double _MIPvalue;/**<MIP Value*/
    // bool _useHistInput;

  public:

    virtual Processor*  newProcessor() { return new ShapingProcessor ; }

    /**Constructor
    */
    ShapingProcessor() ;

    /**Destructor
    */
    ~ShapingProcessor() {};

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

    /** Print Parameters
    */
    virtual void printParameters();
    
    //Fabricio
    struct Imp {
      double t;
      double A;
    };
    static double shaper_old(double* x, double* p);
    static double shaper(double x, double tau, int n, vector<Imp> Imps);
    static double input(double* x, double* p);
    //virtual void PulseShapeCRRC(vector<Imp> Imps, bool plotHit);
    virtual void PulseShapeCRRC(vector<Imp> Imps);
    //virtual void EnergyLinearity();
    vector<struct Imp> histInput(vector<Imp> Imps);

  protected:
    std::string _calorimInpCollection;/**<input collection name*/
    std::string _calorimOutCollection;/**<output collection name*/

    vector<std::string> _calorimInpCollections;/**<input collection name*/
    vector<std::string> _calorimOutCollections;/**<output collection name*/
    vector<std::string> _MIP2GeVFactors;
    vector<float> _MIP2GeVFactors_float;
    float _MIP2GeV;
    float _deltaZ;
    vector<std::string> _FixedPosZ;
    vector<float> _FixedPosZ_float;
    int _NSlabs;
    float _FirstSlabPosZ;
    float _SlabSpacing;
    
    // std::string _encoding;

    // float energyThreshold;/**<Threshold*/
    // float halfdS;/**<half dead space*/

    int _nRun ;/**<run number*/
    int _nEvt ;/**<event number*/
    int _nLeft;/**<number of hits filtered out by shaping*/
    int _nTot;/**<number of total hits */
    bool _plotHit, _plotted;
    //int _nThr;
    //vector<double> _Thrs, _FSTime, _SSRead;
    double _MIPThreshold, _threshold, _FSTime, _SSRead, _FSNoise, _SSNoise;

    double _bwI, _tauF, _tauS, _delay, _bwF, _bwS;
    int _nbinsI, _nbinsS, _nbinsF, _nF, _nS;
    bool _useHistInput, _filterNoise;
    std::string _shapingAuxFilename;
    TFile* _rootout;
    TTree* _treeout;
    double energy_sha, time_sha, energy_sim, time_sim, t_maxFS, t_maxSS, maxFS, maxSS;
    double posx, posy, posz;

  } ;
}
#endif
