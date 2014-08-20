#ifndef CATTools_Muon_H
#define CATTools_Muon_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// Define typedefs for convenience
namespace cat {
  class Muon;
  typedef std::vector<Muon>              MuonCollection;
  typedef edm::Ref<MuonCollection>       MuonRef;
  typedef edm::RefVector<MuonCollection> MuonRefVector;
}

namespace cat {

  class Muon : public Particle{
  public:
    Muon();
    Muon(const reco::LeafCandidate & aMuon); 
    virtual ~Muon();

    double chargedHadronIso(double dR=0.3) const { 
      if( dR == 0.3){ return chargedHadronIso03_;}
      else if( dR == 0.4) { return chargedHadronIso04_;}
      else return -1.0;  
    }
    double puChargedHadronIso(double dR=0.3) const { 
      if( dR == 0.3) { return puChargedHadronIso03_;}
      else if( dR == 0.4) { return puChargedHadronIso04_;}
      else return -1.0;
    }
    double neutralHadronIso(double dR=0.3) const { 
      if( dR == 0.3) { return neutralHadronIso03_;}
      else if( dR == 0.4) { return neutralHadronIso04_;}
      else return -1.0;
    }
    double photonIso(double dR=0.3) const { 
      if( dR == 0.3) { return photonIso03_;}
      else if( dR == 0.4) { return photonIso04_;} 
      else return -1.0;
    }

    double absIso(double dR=0.3, double dBetaFactor=0.5) const{

      if(dBetaFactor>0 && puChargedHadronIso(dR)<0) return -1;

      double neutralIso = neutralHadronIso(dR) + photonIso(dR);
      double corNeutralIso = neutralIso - dBetaFactor * puChargedHadronIso(dR);

      double charged = chargedHadronIso(dR);

      return charged + ( corNeutralIso>0 ? corNeutralIso : 0 ) ;
    }

    double relIso(double dR=0.3, double dBetaFactor=0.5) const{
      double abs = absIso(dR, dBetaFactor)/this->pt();
      return abs >=0 ? abs : -1;
    }

    bool isGlobalMuon() const { return isGlobalMuon_; }
    bool isPFMuon() const { return isPFMuon_; }
    bool isTightMuon() const { return isTightMuon_; }
    bool isLooseMuon() const { return isLooseMuon_; } 
    bool isSoftMuon() const { return isSoftMuon_; } 

    double normalizedChi2() { return normalizedChi2_; }
    int numberOfValidHits() { return numberOfValidHits_; }
    int numberOfValidMuonHits() { return numberOfValidMuonHits_; }
    int numberOfMatchedStations() { return numberOfMatchedStations_; } 
    int numberOfValidPixelHits() { return numberOfValidPixelHits_; }
    int trackerLayersWithMeasurement() { return trackerLayersWithMeasurement_; }
 
    double dxy(){ return dxy_; }
    double dz(){ return dz_; }

    void setChargedHadronIso03(double i) { chargedHadronIso03_ = i; }
    void setPUChargedHadronIso03(double i) { puChargedHadronIso03_ = i; }
    void setNeutralHadronIso03(double i) { neutralHadronIso03_ = i; }
    void setPhotonIso03(double i) { photonIso03_ = i; }

    void setChargedHadronIso04(double i) { chargedHadronIso04_ = i; }
    void setPUChargedHadronIso04(double i) { puChargedHadronIso04_ = i; }
    void setNeutralHadronIso04(double i) { neutralHadronIso04_ = i; }
    void setPhotonIso04(double i) { photonIso04_ = i; }

    void setIsGlobalMuon(bool d) { isGlobalMuon_ = d; }
    void setIsPFMuon(bool d) { isPFMuon_ = d; }
    void setIsTightMuon(bool d) { isTightMuon_ = d; }
    void setIsLooseMuon(bool d) { isLooseMuon_ = d; }
    void setIsSoftMuon(bool d) { isSoftMuon_ = d; }
    
    void setNormalizedChi2(double d) { normalizedChi2_ = d; }
    void setNumberOfValidHits(int i) { numberOfValidHits_ = i; }
    void setNumberOfValidMuonHits(int i) { numberOfValidMuonHits_ = i; }
    void setNumberOfMatchedStations(int i) { numberOfMatchedStations_ = i; }
    void setNumberOfValidPixelHits(int i) { numberOfValidPixelHits_ = i; }
    void setTackerLayersWithMeasurement(int i) { trackerLayersWithMeasurement_ = i; }
 
    void setDz(int d) { dz_ = d; }
    void setDxy(int d) { dxy_ = d; }

  private:

    double chargedHadronIso03_;
    double puChargedHadronIso03_;
    double neutralHadronIso03_;
    double photonIso03_;

    double chargedHadronIso04_;
    double puChargedHadronIso04_;
    double neutralHadronIso04_;
    double photonIso04_;

    bool isGlobalMuon_; 
    bool isPFMuon_; 
    bool isTightMuon_; 
    bool isLooseMuon_; 
    bool isSoftMuon_;

    double normalizedChi2_; 
    int numberOfValidHits_; 
    int numberOfValidMuonHits_; 
    int numberOfMatchedStations_; 
    int numberOfValidPixelHits_; 
    int trackerLayersWithMeasurement_; 

    double dz_;
    double dxy_;

  };
}

#endif
