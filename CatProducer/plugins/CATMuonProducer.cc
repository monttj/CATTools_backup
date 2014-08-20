/**
  \class    cat::CATMuonProducer CATMuonProducer.h "CATTools/CatProducer/interface/CATMuonProducer.h"
  \brief    CAT Muon 
*/


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FWCore/Utilities/interface/isFinite.h"

namespace cat {

  class CATMuonProducer : public edm::EDProducer {
    public:
      explicit CATMuonProducer(const edm::ParameterSet & iConfig);
      virtual ~CATMuonProducer() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      edm::EDGetTokenT<edm::View<pat::Muon> > src_;
      edm::EDGetTokenT<edm::View<reco::Vertex> > vertexLabel_;
      edm::EDGetTokenT<reco::BeamSpot> beamLineSrc_;

  };

} // namespace

cat::CATMuonProducer::CATMuonProducer(const edm::ParameterSet & iConfig) :
    src_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("src"))),
    vertexLabel_(consumes<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
    beamLineSrc_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamLineSrc")))
{
    produces<std::vector<cat::Muon> >();
}

void 
cat::CATMuonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<pat::Muon> > src;
    iEvent.getByToken(src_, src);

    Handle<View<reco::Vertex> > recVtxs;
    iEvent.getByToken(vertexLabel_,recVtxs);

    Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamLineSrc_, beamSpotHandle);

    reco::Vertex pv = recVtxs->at(0);
   
    reco::BeamSpot beamSpot = *beamSpotHandle;
    reco::TrackBase::Point beamPoint(0,0,0);
    beamPoint = reco::TrackBase::Point ( beamSpot.x0(), beamSpot.y0(), beamSpot.z0() );  
 
    auto_ptr<vector<cat::Muon> >  out(new vector<cat::Muon>());

    for (View<pat::Muon>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
      unsigned int idx = it - src->begin();
      const pat::Muon & aPatMuon = src->at(idx);

      cat::Muon aMuon(aPatMuon);

      aMuon.setChargedHadronIso04( aPatMuon.chargedHadronIso() );
      aMuon.setNeutralHadronIso04( aPatMuon.neutralHadronIso() );
      aMuon.setPhotonIso04( aPatMuon.photonIso() );
      aMuon.setPUChargedHadronIso04( aPatMuon.puChargedHadronIso() );

      aMuon.setChargedHadronIso03( aPatMuon.userIsolation("pat::User1Iso") );
      aMuon.setNeutralHadronIso03( aPatMuon.userIsolation("pat::User2Iso") );
      aMuon.setPhotonIso03( aPatMuon.userIsolation("pat::User3Iso") );
      aMuon.setPUChargedHadronIso03( aPatMuon.userIsolation("pat::User4Iso") );

      aMuon.setIsGlobalMuon( aPatMuon.isGlobalMuon() );
      aMuon.setIsPFMuon( aPatMuon.isPFMuon() );
      aMuon.setIsTightMuon( aPatMuon.isTightMuon(pv) );
      aMuon.setIsLooseMuon( aPatMuon.isLooseMuon() );
      aMuon.setIsSoftMuon( aPatMuon.isSoftMuon(pv) );

      aMuon.setNumberOfMatchedStations( aPatMuon.numberOfMatchedStations() );

      if ( aPatMuon.globalTrack().isNonnull() && aPatMuon.globalTrack().isAvailable() ) {
        aMuon.setNormalizedChi2( aPatMuon.normChi2() );
        aMuon.setNumberOfValidMuonHits( aPatMuon.globalTrack()->hitPattern().numberOfValidMuonHits() );
      }

      if ( aPatMuon.innerTrack().isNonnull() && aPatMuon.innerTrack().isAvailable() ){
        aMuon.setNumberOfValidHits( aPatMuon.numberOfValidHits() );
        aMuon.setNumberOfValidPixelHits( aPatMuon.innerTrack()->hitPattern().numberOfValidPixelHits() );
        aMuon.setTackerLayersWithMeasurement( aPatMuon.innerTrack()->hitPattern().trackerLayersWithMeasurement() ); 
      }

      double dxy = fabs(aPatMuon.muonBestTrack()->dxy(pv.position()));
      aMuon.setDxy( dxy );
      double dz = fabs(aPatMuon.muonBestTrack()->dz(pv.position()));
      aMuon.setDz( dz ); 
 
      out->push_back(aMuon);

    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATMuonProducer);
