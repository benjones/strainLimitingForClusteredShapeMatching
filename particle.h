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
  double totalweight; // sum of weights of this particle in all clusters
  size_t numClusters; // only used to determine singelton particles
  bool outsideSomeMovingPlane;
  std::vector<int> clusters; 
  Particle() {};
  Particle(const Particle &p)  = default; //c++11 magic :)
  //: position(p.position), velocity(p.velocity), restPosition(p.restPosition),
  //	oldPosition(p.oldPosition), goalPosition(p.goalPosition), goalVelocity(p.goalVelocity), mass(p.mass),
  //	totalweight(p.totalweight), numClusters(p.numClusters), outsideSomeMovingPlane(p.outsideSomeMovingPlane),
  //	clusters(p.clusters){};
};

class Cluster {
 public:
  Eigen::Vector3d restCom;
  Eigen::Vector3d worldCom;
  Eigen::Matrix3d aInv, Fp, FpNew;
  std::vector<std::pair<int, double> > neighbors;
  double mass, width, renderWidth;
  double toughness;
  double cstrain; // cumulative strain (for work hardening)
  Cluster() {Fp.setIdentity(); cstrain = 0.0;}
  void updatePlasticity(const Eigen::Vector3d &sigma, const Eigen::Matrix3d &U, const Eigen::Matrix3d &V, 
	  double yield, double nu, double hardening);
};
