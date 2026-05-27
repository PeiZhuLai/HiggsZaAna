#ifndef RecoEgamma_EgammaMLPhotonProducers_MLPhotonProducer_h
#define RecoEgamma_EgammaMLPhotonProducers_MLPhotonProducer_h

#include <iostream>
#include <vector>
#include <memory>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
// #include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "MLDataFormats/EgammaCandidates/interface/MLPhoton.h"
#include "MLDataFormats/EgammaCandidates/interface/MLPhotonFwd.h"

#include "RecoEgamma/EgammaMLPhotonProducers/interface/Cluster.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"


// #include "Geometry/CaloGeometry/interface/CaloGeometry.h"
// #include "Geometry/CaloTopology/interface/CaloTopology.h"
// #include "DataFormats/EgammaCandidates/interface/Conversion.h"
// #include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
// #include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
// #include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
// #include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
// #include "Geometry/Records/interface/CaloGeometryRecord.h"
// #include "Geometry/CaloGeometry/interface/CaloGenericDetId.h"
// #include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"

// using namespace cms::Ort;
// using namespace edm;
// using namespace fastjet;
// using namespace fastjet::contrib;

class MLPhotonProducer : public edm::stream::EDProducer<> {
   public:
      MLPhotonProducer( const edm::ParameterSet& ps);
      ~MLPhotonProducer() override;

      void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;


   private:
      // ----------member functions ---------------------------
      math::PtEtaPhiMLorentzVector calculateLorentzVector(float moe, float energy, float eta, float phi, float zpv);
      std::pair<int, float> FindNearest(Cluster C1, std::vector<Cluster> inC);
      static bool comparator(const Cluster& C1, const Cluster& C2);
      std::vector<Cluster> DoPairings(std::vector<Cluster> inC, float R);

      // ----------member data ---------------------------
      std::string collectionLabel;

      // ML
      cms::Ort::ONNXRuntime ortClassifier;
      cms::Ort::ONNXRuntime ortRegressor;

      const double MATCH_DR = 0.15; //Cluster Size
      const double ISO_DR = 0.3; //Isolation Cone Size
      edm::EDGetTokenT<std::vector<reco::CaloCluster>> clustersToken;
      edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> HEEToken;
      edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> HEBToken;

      const edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken;
      const edm::EDGetTokenT<std::vector<pat::PackedCandidate>> pfCandToken;


  };
  #endif