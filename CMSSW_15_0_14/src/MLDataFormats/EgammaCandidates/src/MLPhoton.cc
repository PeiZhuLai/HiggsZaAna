#include "MLDataFormats/EgammaCandidates/interface/MLPhoton.h"

using namespace reco;

MLPhoton::MLPhoton(const Candidate::PolarLorentzVector& p4, const Point vtx):
	RecoCandidate(0, p4, vtx, 0),
	p4_(p4) {}

bool MLPhoton::overlap(const Candidate& c) const {
	return reco::deltaR(p4_, c.p4()) < 0.1;
}

MLPhoton::~MLPhoton() {}