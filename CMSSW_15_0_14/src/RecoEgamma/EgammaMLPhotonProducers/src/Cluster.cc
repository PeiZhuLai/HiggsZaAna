#include "RecoEgamma/EgammaMLPhotonProducers/interface/Cluster.h"

Cluster::Cluster(const float &Eta, const float &Phi, const int &iEta, const int &iPhi, const float &E, const bool &nCrack){
	vec.SetCoordinates(E/std::cosh(Eta), Eta, Phi);
	Es.push_back(E);	
	ietas.push_back(iEta);
	iphis.push_back(iPhi);
	ncracks.push_back(nCrack);
}

Cluster::Cluster(const Cluster & C){
	vec = C.vec;
	Es = C.Es;
	ietas = C.ietas;
	iphis = C.iphis;
	ncracks = C.ncracks;
}

Cluster::~Cluster(void){}

void Cluster::combine(Cluster C){
	vec+=C.vec;
	ietas.insert( ietas.end(), C.ietas.begin(), C.ietas.end() );
	iphis.insert( iphis.end(), C.iphis.begin(), C.iphis.end() );
	ncracks.insert( ncracks.end(), C.ncracks.begin(), C.ncracks.end() );
	Es.insert( Es.end(), C.Es.begin(), C.Es.end() );
}

void Cluster::makeImage(){
	int iphi_max = *std::max_element(iphis.begin(), iphis.end());
	int iphi_min = *std::min_element(iphis.begin(), iphis.end());

	if(iphi_max - iphi_min > isize && iphi_max > 340 && iphi_min < 20){  //This loop fixes wrapped events
	for(unsigned int ii=0; ii<iphis.size(); ii++){
		if(iphis[ii] < isize){
		iphis[ii] = 360 + iphis[ii] + 1;
		}
	}
	}

	//Add minimum + half isize to each element of ieta, iphi
	int ieta_min = *std::min_element(ietas.begin(), ietas.end());
	iphi_min = *std::min_element(iphis.begin(), iphis.end());

	for(int& d : ietas){d += ieta_min + (isize/2) ;}
	for(int& d : iphis){d += iphi_min + (isize/2) ;}

	int cxval = ietas[0];
	int cyval = iphis[0];

	//std::vector<int> xcoords;
	//std::vector<int> ycoords;

	for(unsigned int ii=0; ii<ietas.size(); ii++){
	xcoords.push_back(ietas[ii] - cxval + (isize / 2));
	ycoords.push_back(iphis[ii] - cyval + (isize / 2));
	}

	int y_min = *std::min_element(ycoords.begin(), ycoords.end());
	while(y_min < 0){
	for(int& d : ycoords){d += abs(y_min);}
	y_min = *std::min_element(ycoords.begin(), ycoords.end());
	}

	int y_max = *std::max_element(ycoords.begin(), ycoords.end());
	if(y_max > isize){
	for(int& d : ycoords){d -= y_max - isize;}
	}

	///////////////////////////
	std::vector<std::vector<float>> img(isize, std::vector<float>(isize)); //initialize vector of size isize x isize to all zeros
	for(unsigned int ii=0; ii<xcoords.size(); ii++){
	int x = xcoords[ii];
	int y = ycoords[ii];
	float e = Es[ii];
	if(x >= isize){x=isize-1;}
	if(y >= isize){y=isize-1;}
	if(x < 0){x=0;}
	if(y < 0){y=0;}
	img[x][y] += e;
	}

	////////////////////////
	//Flatten Vector
	///////////////////////

	std::vector<float>flatImg;//new 1D vector
	// using nested for loops to construct a new 1D vector
	for(unsigned int i=0; i<img.size(); i++){
		for(unsigned int j=0; j<img[i].size(); j++){
			flatImg.push_back(img[i][j]);
		}
	}

	for(unsigned int i=0; i<isize*isize; i++){
		image.push_back(flatImg);
	}

	return;
}

float Cluster::compute_En(float n){
	float totalE = 0.;
	float numerator = 0.;
	totalE = std::accumulate(Es.begin(), Es.end(), 0.0);
	for(unsigned int ii=0; ii < xcoords.size(); ii++){
		int xx = xcoords.at(ii);
		int yy = ycoords.at(ii);
		numerator += pow( (xx * xx + yy * yy), n/2 ) ;
	}

	return numerator / totalE;
}