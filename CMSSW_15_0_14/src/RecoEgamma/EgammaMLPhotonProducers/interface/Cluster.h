
#include <numeric>
#include <algorithm>

#include <cmath>
#include <vector>

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/deltaR.h"

class Cluster {
	public:

    	static const int isize=30;

		math::RhoEtaPhiVectorF vec;
		std::vector<int> ietas;
		std::vector<int> iphis;
		std::vector<int> ncracks;
		std::vector<float> Es;
		std::vector<int> xcoords;
		std::vector<int> ycoords;
		std::vector<std::vector<float>> image;
		//std::vector<std::vector<float>> image(isize*isize, std::vector<float>(isize*isize,0.0));
		Cluster(const float &Eta, const float &Phi, const int &iEta, const int &iPhi, const float &E, const bool &nCrack);
		Cluster(const Cluster & clu);

		~Cluster();


		float Eta(){return vec.Eta();}
		float Phi(){return vec.Phi();}

		float deltaR(Cluster C){return reco::deltaR(vec, C.vec);}

		void combine(Cluster C);
		float getTotalE(void){return std::accumulate(Es.begin(), Es.end(), 0.0f);}
		void makeImage();
		float compute_En(float);

};