#ifndef EgammaReco_MLPhotonFwd_h
#define EgammaReco_MLPhotonFwd_h

#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"

namespace reco {
  class MLPhoton;

  /// collectin of MLPhoton objects
  typedef std::vector<MLPhoton> MLPhotonCollection;

  /// reference to an object in a collection of MLPhoton objects
  typedef edm::Ref<MLPhotonCollection> MLPhotonRef;

  /// reference to a collection of MLPhoton objects
  typedef edm::RefProd<MLPhotonCollection> MLPhotonRefProd;

  /// vector of objects in the same collection of MLPhoton objects
  typedef edm::RefVector<MLPhotonCollection> MLPhotonRefVector;

  /// iterator over a vector of reference to MLPhoton objects
  typedef MLPhotonRefVector::iterator ml_photon_iterator;
}  // namespace reco

#endif
