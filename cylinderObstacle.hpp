#pragma once

#include <Eigen/Dense>

class Particle;

class CylinderObstacle{

public:
  Eigen::Vector3d normal, supportPoint;
  double radius;
  CylinderObstacle(const Eigen::Vector3d& _normal, const Eigen::Vector3d& _supportPoint, double _radius)
	:normal(_normal.normalized()), supportPoint{_supportPoint}, radius{_radius}
  {}

  void bounceParticle(Particle& particle) const;
};
