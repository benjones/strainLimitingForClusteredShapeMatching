#pragma once
#include <Eigen/Dense>
#include <vector>

class Particle{
public:
  Eigen::Vector3d position, 
	velocity, 
	restPosition, 
	oldPosition, 
	goalPosition, 
	goalVelocity;
  double mass; 
  size_t numClusters;
  bool outsideSomeMovingPlane;
  std::vector<int> clusters;
  std::vector<double> weights;
};

class Cluster {
 public:
	Eigen::Vector3d restCom;
  Eigen::Vector3d worldCom;
  Eigen::Matrix3d aInv, Fp, FpNew;
	std::vector<int> neighbors;
  double mass, width, renderWidth;
  double toughness;
  double cstrain; // cumulative strain (for work hardening)
  Cluster() {Fp.setIdentity(); cstrain = 0.0;}
};
