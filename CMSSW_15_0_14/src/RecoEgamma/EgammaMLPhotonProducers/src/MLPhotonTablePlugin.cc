// Register a SimpleFlatTableProducer template specialised for reco::MLPhoton
// so the NanoAOD machinery can expand member-function expressions like
// "diphotonScore()" / "massEnergyRatio()" defined on the derived class.

#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"
#include "MLDataFormats/EgammaCandidates/interface/MLPhoton.h"

typedef SimpleFlatTableProducer<reco::MLPhoton> SimpleMLPhotonFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SimpleMLPhotonFlatTableProducer);
