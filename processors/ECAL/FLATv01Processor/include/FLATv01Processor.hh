#ifndef FLATv01Processor_HH
#define FLATv01Processor_HH 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <map>

#include <EVENT/MCParticle.h>

//ROOT
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

using namespace lcio ;
using namespace marlin ;
using std::vector;

namespace CALICE
{
  /**
  * @brief Class doing the conversion between GeV and MIP + Correction of K layer (not smart!!)
  *
  * @autho fabricio.jm@cern.ch, based on that of
  * eldwan.brianne@desy.de
  * @version 1.0
  */

  class FLATv01Processor : public Processor {

  private:
    //processor parameters:
    // float _MIPvalue; /**<MIP conversion factor in GeV*/
    std::string _rootfilename; /**<rootfile name*/
    std::string _eConfName; /**<econf name*/

  public:

    virtual Processor*  newProcessor() { return new FLATv01Processor ; }

    /**Constructor
    */
    FLATv01Processor() ;

    /**Destructor
    */
    ~FLATv01Processor() {};

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


  protected:
    std::string _calorimInpCollection; /**<input collection name*/
    std::string _calorimOutCollection; /**<output collection name*/

    vector<std::string> _calorimInpCollections;/**<input collection name*/
    vector<std::string> _calorimOutCollections;/**<output collection name*/

    std::map<std::string, int> _colDict;

    int _nRun ; /**<run number*/
    int _nEvt ; /**<evt number*/
    // Event
    int _nHits;
    vector<MCParticle *> _particles;
    //std::vector<int> _nSubHits;
    std::vector<float> _energy, _positionX, _positionY, _positionZ;
    std::vector<int> _cellID0, _cellID1, _nMCContributions, _colType;
    std::vector<float> _maxE_subhit_energy, _maxE_subhit_time;
    std::vector<float> _maxE_subhit_positionX, _maxE_subhit_positionY, _maxE_subhit_positionZ;
    std::vector<float> _maxE_particle_pX, _maxE_particle_pY, _maxE_particle_pZ;
    std::vector<float> _maxE_particle_energy;
    std::vector<int> _maxE_subhit_PDG, _maxE_particle_PDG, _maxE_particle_index;
    //std::vector<std::vector<float>> _energyCont, _timeCont, _stepPositionX, _stepPositionY, _stepPositionZ;
    //std::vector<std::vector<float>> _PDGCont;
    
    TFile* _rootout;
    TTree* _treeout;

  } ;
}
#endif
