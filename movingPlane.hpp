#pragma once
#include <Eigen/Dense>
class Particle;

class MovingPlane{
public:

  MovingPlane(Eigen::Vector3d _normal, double _offset, double _velocity)
	:normal(_normal.normalized()), offset(_offset), velocity(_velocity)
  {}

  void bounceParticle(Particle& particle, double timeElapsed) const;

  Eigen::Vector3d normal;
  double offset, velocity;
};
