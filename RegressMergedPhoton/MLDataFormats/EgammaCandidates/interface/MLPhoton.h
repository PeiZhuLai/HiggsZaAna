// Description: Class to store ML photon information
//

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

// #include "DataFormats/Math/interface/LorentzVector.h"

namespace reco {

	class MLPhoton : public RecoCandidate {
		
		public:

			// Default constructor
			MLPhoton() : RecoCandidate() {}

			MLPhoton(const PolarLorentzVector& p4, const Point vtx);

			~MLPhoton();

			// Getter methods
			float massEnergyRatio() const {return moe_;}

			float diphotonScore() const {return diphoton_score_;}
			float monophotonScore() const {return monophoton_score_;}
			float hadronScore() const {return hadron_score_;}

			float pfIsolation() const {return pfIso_;}

			float r1() const {return r1_;}
			float r2() const {return r2_;}
			float r3() const {return r3_;}

			// Setter methods
			void setMassEnergyRatio(float moe) {moe_ = moe;}

			void setDiphotonScore(float score) {diphoton_score_ = score;}
			void setMonophotonScore(float score) {monophoton_score_ = score;}
			void setHadronScore(float score) {hadron_score_ = score;}

			void setPFIsolation(float iso) {pfIso_ = iso;}

			void setR1(float r1) {r1_ = r1;}
			void setR2(float r2) {r2_ = r2;}
			void setR3(float r3) {r3_ = r3;}


		private:

			bool overlap(const Candidate&) const override;		
			
			PolarLorentzVector p4_;

			float moe_;
			float energy_;
			float eta_;
			float phi_;

			float diphoton_score_;
			float monophoton_score_;
			float hadron_score_;

			float pfIso_;

			float r1_;
			float r2_;
			float r3_;

			// float pfIso;

	};

}