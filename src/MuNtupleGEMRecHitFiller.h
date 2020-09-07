#ifndef MuNtuple_MuNtupleGEMRecHitFiller_h
#define MuNtuple_MuNtupleGEMRecHitFiller_h


#include "MuDPGAnalysis/MuonDPGNtuples/src/MuNtupleBaseFiller.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"

#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"

#include <vector>

class MuNtupleGEMRecHitFiller : public MuNtupleBaseFiller
{
 public:

  //Constructor
  MuNtupleGEMRecHitFiller(edm::ConsumesCollector && collector,
			  const std::shared_ptr<MuNtupleConfig> config,
			  std::shared_ptr<TTree> tree, const std::string & label);


  //Destructor
  virtual ~MuNtupleGEMRecHitFiller();

  virtual void initialize() final; 
  
  virtual void clear() final;

  virtual void fill(const edm::Event & ev) final;

 
 private:

  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> geomToken_;
  edm::EDGetTokenT<GEMRecHitCollection> m_rechit_token_;
  
  unsigned int m_nRecHits;
  
  std::vector<float>  m_rechit_g_r;
  std::vector<float>  m_rechit_g_phi;
  std::vector<float>  m_rechit_g_x;
  std::vector<float>  m_rechit_g_y;
  std::vector<float>  m_rechit_g_z;
  
};

#endif
