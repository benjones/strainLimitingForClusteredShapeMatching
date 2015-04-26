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
  double totalweight; // only used for initialization
  size_t numClusters; // only used to determine singelton particles
  bool outsideSomeMovingPlane;
  std::vector<int> clusters; 
  Particle() {};
  Particle(const Particle &p) : position(p.position), velocity(p.velocity), restPosition(p.restPosition),
	oldPosition(p.oldPosition), goalPosition(p.goalPosition), goalVelocity(p.goalVelocity), mass(p.mass),
	totalweight(p.totalweight), numClusters(p.numClusters), outsideSomeMovingPlane(p.outsideSomeMovingPlane),
	clusters(p.clusters){};
};

class Cluster {
 public:
	Eigen::Vector3d restCom;
  Eigen::Vector3d worldCom;
  Eigen::Matrix3d aInv, Fp, FpNew;
	std::vector<int> neighbors;
	std::vector<double> weights;
  double mass, width, renderWidth;
  double toughness;
  double cstrain; // cumulative strain (for work hardening)
  Cluster() {Fp.setIdentity(); cstrain = 0.0;}
};
