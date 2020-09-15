#include "MuDPGAnalysis/MuonDPGNtuples/src/MuNtupleGEMMuonFiller.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Common/interface/Ref.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegment.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include <iostream>

#include "TVectorF.h"

MuNtupleGEMMuonFiller::MuNtupleGEMMuonFiller(edm::ConsumesCollector && collector,
				       const std::shared_ptr<MuNtupleConfig> config, 
				       std::shared_ptr<TTree> tree, const std::string & label) : 
  MuNtupleBaseFiller(config, tree, label), m_nullVecF()
{

  edm::InputTag & muonTag = m_config->m_inputTags["muonTag"];
  if (muonTag.label() != "none") m_muToken = collector.consumes<reco::MuonCollection>(muonTag);

  edm::InputTag & primaryVerticesTag = m_config->m_inputTags["primaryVerticesTag"];
  if (primaryVerticesTag.label() != "none") m_primaryVerticesToken = collector.consumes<std::vector<reco::Vertex>>(primaryVerticesTag);

  edm::InputTag & gemSegmentTag = m_config->m_inputTags["gemSegmentTag"];
  if (gemSegmentTag.label() != "none") m_gemSegmentToken = collector.consumes<GEMSegmentCollection>(gemSegmentTag);
  
  edm::InputTag & cscSegmentTag = m_config->m_inputTags["cscSegmentTag"];
  if(cscSegmentTag.label() != "none") m_cscSegmentToken = collector.consumes<CSCSegmentCollection>(cscSegmentTag);

  edm::InputTag & trigResultsTag = m_config->m_inputTags["trigResultsTag"];
  if (trigResultsTag.label() != "none") m_trigResultsToken = collector.consumes<edm::TriggerResults>(trigResultsTag);

  edm::InputTag & trigEventTag = m_config->m_inputTags["trigEventTag"];
  if (trigEventTag.label() != "none") m_trigEventToken = collector.consumes<trigger::TriggerEvent>(trigEventTag);

  edm::InputTag & gemRecHitTag = m_config->m_inputTags["gemRecHitTag"];
  if (gemRecHitTag.label() != "none") m_gemRecHitToken = collector.consumes<GEMRecHitCollection>(gemRecHitTag);
  
}

MuNtupleGEMMuonFiller::~MuNtupleGEMMuonFiller() 
{ 

}

void MuNtupleGEMMuonFiller::initialize()
{

  m_tree->Branch((m_label + "_nMuons").c_str(), &m_nMuons);
  
  m_tree->Branch((m_label + "_pt").c_str(), &m_pt);
  m_tree->Branch((m_label + "_phi").c_str(), &m_phi);
  m_tree->Branch((m_label + "_eta").c_str(), &m_eta);
  m_tree->Branch((m_label + "_charge").c_str(), &m_charge);

  m_tree->Branch((m_label + "_isGlobal").c_str(), &m_isGlobal);
  m_tree->Branch((m_label + "_isStandalone").c_str(), &m_isStandalone);
  m_tree->Branch((m_label + "_isTracker").c_str(), &m_isTracker);
  //m_tree->Branch((m_label + "_isTrackerArb").c_str(), &m_isTrackerArb);
  m_tree->Branch((m_label + "_isGEM").c_str(), &m_isGEM);

  m_tree->Branch((m_label + "_isLoose").c_str(), &m_isLoose);
  m_tree->Branch((m_label + "_isMedium").c_str(), &m_isMedium);
  m_tree->Branch((m_label + "_isTight").c_str(), &m_isTight);

}

void MuNtupleGEMMuonFiller::clear()
{

  m_nMuons = 0;

  m_pt.clear();
  m_phi.clear();
  m_eta.clear();
  m_charge.clear();

  m_isGlobal.clear();
  m_isStandalone.clear();
  m_isTracker.clear();
  //m_isTrackerArb.clear();
  m_isGEM.clear();

  m_isLoose.clear();
  m_isMedium.clear();
  m_isTight.clear();

}

void MuNtupleGEMMuonFiller::fill(const edm::Event & ev)
{

  clear();

  auto muons = conditionalGet<reco::MuonCollection>(ev, m_muToken, "MuonCollection");
  auto gem_segments = conditionalGet<GEMSegmentCollection>(ev,m_gemSegmentToken, "GEMSegmentCollection");
  auto csc_segments = conditionalGet<CSCSegmentCollection>(ev,m_cscSegmentToken, "CSCSegmentCollection" );
  auto vtxs = conditionalGet<std::vector<reco::Vertex>>(ev, m_primaryVerticesToken, "std::vector<reco::Vertex>");

  edm::Handle<GEMRecHitCollection> rechit_collection;
  ev.getByToken(m_gemRecHitToken, rechit_collection);
  if (not rechit_collection.isValid()) {
    std::cout << "GEMRecHitCollection is invalid" << std::endl;
    return;
  }
  
  edm::ESHandle<Propagator>&& propagator = m_config->muon_service->propagator("SteppingHelixPropagatorAny");
  if (not propagator.isValid()) {
    std::cout<< "Propagator is invalid" << std::endl;
    return;
  }

  
  const auto gem = m_config->m_gemGeometry;
  if (not gem.isValid()) {
    std::cout << "GEMGeometry is invalid" << std::endl;
    return;
  }
  
  const auto transient_track_builder = m_config->m_transientTrackBuilder;


  if (muons.isValid() && csc_segments.isValid() && vtxs.isValid()) 
    {
      //std::cout << "ciao" <<std::endl;
      for (const auto & muon : (*muons))
	{

	  m_pt.push_back(muon.pt());
	  m_eta.push_back(muon.eta());
	  m_phi.push_back(muon.phi());
	  m_charge.push_back(muon.charge());
	  
	  m_isGlobal.push_back(muon.isGlobalMuon());
	  m_isStandalone.push_back(muon.isStandAloneMuon());
	  m_isTracker.push_back(muon.isTrackerMuon());
	  //m_isTrackerArb.push_back(isTrackerArb());
	  //m_isGEM.push_back(muon.isGEM());

	  m_isLoose.push_back(muon.passed(reco::Muon::CutBasedIdLoose));
	  m_isMedium.push_back(muon.passed(reco::Muon::CutBasedIdMedium));
	  m_isTight.push_back(muon.passed(reco::Muon::CutBasedIdTight));

	  m_nMuons++;

	  if(!muon.outerTrack().isNull())
	    {
	      //reco::TrackRef outerTrackRef = muon.outerTrack();
	      //std::cout << "ciao" <<std::endl;
	      //auto recHitMu = outerTrackRef->recHitsBegin();
	      //auto recHitMuEnd = outerTrackRef->recHitsEnd();

	      /*for(; recHitMu != recHitMuEnd; ++recHitMu)
		{
		  DetId detId = (*recHitMu)->geographicalId();
		  //std::cout << "ciao" <<std::endl;
		  if(detId.det() == DetId::Muon && detId.subdetId() == MuonSubdetId::CSC)
		    {
		      //std::cout << (*recHitMu)->chamberId().station() <<std::endl;
		      const auto cscSegmentSta = dynamic_cast<const CSCSegment*>((*recHitMu));
		      //std::cout << (*recHitMu)->localPosition().x() << std::endl;
		      //std::cout << "!!" << std::endl;
		      //std::cout << "ciao" <<std::endl;

		      for(const auto& csc_segment : (*csc_segments))
			{
			  //std::cout << csc_segment.localPosition().x() <<std::endl;
			  //std::cout << "ciao" <<std::endl;
			  if( 
			     std::abs((*recHitMu)->localPosition().x() -  csc_segment.localPosition().x()) < 0.01
			     //cscSegmentSta &&
			    //cscSegmentSta->cscDetId().station() == csc_segment.cscDetId().station()
			  )
			  {
			    //std::cout << "ciao" << std::endl; 
			   }
		       
			}
		    }      
		  
		}*/
	      
	      const reco::Track* track = nullptr;
	      track = muon.outerTrack().get();
	      if (track == nullptr) {
		std::cout << "failed to get muon track" << std::endl;
		continue;
	      }
	      
	      const reco::TransientTrack&& transient_track = transient_track_builder->build(track);
	      if (not transient_track.isValid()) {
		std::cout<<"failed  to build TransientTrack" << std::endl;
		continue;
	      }

	      for (const GEMEtaPartition* eta_partition : gem->etaPartitions()) {
		// Skip propagation inn the opposite direction.
		if (muon.eta() * eta_partition->id().region() < 0)
		  continue;

		const BoundPlane& bound_plane = eta_partition->surface();
		
		const TrajectoryStateOnSurface&& tsos =
		  propagator->propagate(transient_track.outermostMeasurementState(), bound_plane);
		if (not tsos.isValid()) {
		  continue;
		}
		
		const LocalPoint&& tsos_local_pos = tsos.localPosition();
		const LocalPoint tsos_local_pos_2d(tsos_local_pos.x(), tsos_local_pos.y(), 0.0f);
		if (not bound_plane.bounds().inside(tsos_local_pos_2d)) {
		  continue;
		}
		
		const GEMDetId&& gem_id = eta_partition->id();
		/*const GEMRecHit* matched_hit = findMatchedHit(tsos_local_pos.x(), rechit_collection->get(gem_id));
		const LocalPoint&& hit_local_pos = matched_hit->localPosition();
		const GlobalPoint&& hit_global_pos = eta_partition->toGlobal(hit_local_pos);
		const GlobalPoint&& tsos_global_pos = tsos.globalPosition();
		*/
		
	      }

	    }
	}
    }

  return;

}

