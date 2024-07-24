#ifndef LCIO2BuildProcessor_HH
#define LCIO2BuildProcessor_HH 1

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

  class LCIO2BuildProcessor : public Processor {

  private:
    //processor parameters:
    // float _MIPvalue; /**<MIP conversion factor in GeV*/
    std::string _rootfilename; /**<rootfile name*/
    std::string _outFileName; /**<econf name*/

  public:

    virtual Processor*  newProcessor() { return new LCIO2BuildProcessor ; }

    /**Constructor
    */
    LCIO2BuildProcessor() ;

    /**Destructor
    */
    ~LCIO2BuildProcessor() {};

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
    
    // static void hitCast();
    // static SimCalorimeterHit hitCast(SimCalorimeterHit* aHit);
    // static SimCalorimeterHit* hitCast(SimCalorimeterHit* aHit);
    //virtual CalorimeterHit* hitCast(CalorimeterHit* aHit);



  protected:
    
    
    std::string _calorimInpCollection; /**<input collection name*/
    //    std::string _calorimOutCollection; /**<output collection name*/

    vector<std::string> _calorimInpCollections;/**<input collection name*/
    //    vector<std::string> _calorimOutCollections;/**<output collection name*/

    //    std::map<std::string, int> _colDict;
    vector<std::string> _siThicknesses;
    vector<float> _siThicknesses_float;
    // std::string _hitType;

    int _nRun ; /**<run number*/
    int _nEvt ; /**<evt number*/
    // Event
    int _nHits;
    //std::vector<float> _energy, _positionX, _positionY, _positionZ;
    int event, spill, cycle, bcid, bcid_first_sca_full, bcid_merge_end, id_run, id_dat, nhit_slab, nhit_chip, nhit_chan, nhit_len;
    float sum_energy, sum_energy_w, sum_energy_lg;
    std::vector<int> hit_slab, hit_chip, hit_chan, hit_sca, hit_adc_high, hit_adc_low, hit_n_scas_filled, hit_isHit, hit_isMasked, hit_isCommissioned, hit_positron, hit_nMC;
    std::vector<float> hit_energy, hit_energy_w, hit_energy_lg, hit_x, hit_y, hit_z;
    //std::vector<int> _cellID0, _cellID1, _nMCContributions, _colType;
    
    std::vector<std::string> _FixedPosZ;
    std::vector<std::string> _GeV2MIPFactors;
    float _GeV2MIP;
    std::vector<std::string> _MapFilenames;
    int _NSlabs;
    float _FirstSlabPosZ, _SlabSpacing, _HalfCenterGap;
    std::vector<int> _SlabMapIndices;
    std::vector<float> _FixedPosZ_float;
    std::vector<float> _GeV2MIPFactors_float ;
    std::vector<std::vector<std::vector<float>>> _maps;
    // std::vector<std::vector<std::vector<float>>> _maps(_MapFilenames.size(), std::vector<std::vector<float>>(1024, std::vector<float>(2)));
    float _deltaZ;
    bool _printType;
    bool _ConversionGeV2MIP;
    TFile* _rootout;
    TTree* _treeout;

  } ;
}
#endif
