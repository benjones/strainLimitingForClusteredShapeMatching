#include "particle.h"
#include "utils.h"
#include <iostream>

inline double sqr (double x) {return x*x;}

void Cluster::updatePlasticity(const Eigen::Vector3d &sigma, const Eigen::Matrix3d &U, const Eigen::Matrix3d &V, 
	double yield, double nu, double hardening) {
  FpNew = Fp;
  if (nu > 0.0) {
	if (sigma(2) >= 1e-4) { // adam says: the second clause is a quick hack to avoid plasticity when sigma is degenerate
	  Eigen::Vector3d FpHat = sigma;
	  //std::cout<<FpHat(0)<<" "<<FpHat(1)<<" "<<FpHat(2)<<" => ";
	  FpHat *= 1.0/cbrt(FpHat(0) * FpHat(1) * FpHat(2));
	  //std::cout<<FpHat(0)<<" "<<FpHat(1)<<" "<<FpHat(2)<<std::endl;
	  double norm = sqrt(sqr(FpHat(0)-1.0) + sqr(FpHat(1)-1.0) + sqr(FpHat(2)-1.0));
	  double local_yield = yield + hardening * cstrain;
	  if (norm > local_yield) {	
		double gamma = std::min(1.0, nu * (norm - local_yield) / norm);
		FpHat(0) = pow(FpHat(0), gamma);
		FpHat(1) = pow(FpHat(1), gamma);
		FpHat(2) = pow(FpHat(2), gamma);
		// update cluster.Fp
		FpNew = FpHat.asDiagonal() * V.transpose() * Fp * V.determinant();
	  } 
	}
	cstrain += sqrt(sqr(sigma(0)-1.0) + sqr(sigma(1)-1.0) + sqr(sigma(2)-1.0));
  }
}

CollisionGeometry & CollisionGeometry::operator= (const CollisionGeometry &that) {
  c = that.c; 
  r = that.r; 
  planes = that.planes;
  //for (auto &p : that.planes) 
  //planes.push_back(std::pair<Eigen::Vector3d, double>(p.first, p.second));
  return *this;
};

Eigen::Matrix4d Cluster::getVisTransform() const {
    Eigen::Matrix4d total_t = Eigen::Matrix4d::Identity();
    //step 3, translate from rest space origin to world 
    total_t.block<3,1>(0,3) << worldCom;
    
    //step 2, rotate about rest-to-world transform
    Eigen::Matrix4d rot = Eigen::Matrix4d::Identity();
    //push the rotation onto a 4x4 glMatrix
    auto polar = utils::polarDecomp(restToWorldTransform);
    rot.block<3,3>(0,0) << polar.first;
    total_t = total_t * rot;

    //step 1, translate rest com to origin
    Eigen::Matrix4d local_t = Eigen::Matrix4d::Identity();
    local_t.block<3,1>(0,3) << -restCom;
    total_t = total_t * local_t;

    return total_t;
}
