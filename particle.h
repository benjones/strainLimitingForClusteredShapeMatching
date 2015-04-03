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
  std::vector<int> clusters;
};

class Cluster {
 public:
	Eigen::Vector3d restCom;
  Eigen::Vector3d worldCom;
	Eigen::Matrix3d aInv;
	std::vector<int> neighbors;
	double mass, width;
};
