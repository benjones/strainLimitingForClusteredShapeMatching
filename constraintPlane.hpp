#pragma once

#include <Eigen/Dense>
#include <vector>

class Particle;

class ConstraintPlane{

public:
  ConstraintPlane(double xValue, const std::vector<Particle>& paticles);
  
  void constrainParticles(std::vector<Particle>& particles) const;
	  
  std::vector<std::pair<int, Eigen::Vector3d> > constraints;
};
