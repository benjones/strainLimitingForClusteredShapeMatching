#pragma once
#include <Eigen/Dense>
class Particle;

class MovingPlane{
public:

  MovingPlane(Eigen::Vector3d _normal, double _offset, double _velocity)
	:normal(_normal.normalized()), offset(_offset), velocity(_velocity)
  {}

  bool outside(Particle& particle) const;
  void bounceParticle(Particle& particle, double timeElapsed) const;
  void dragParticle(Particle& particle, double timeElapsed) const;
  void backsideReflectBounceParticle(Particle& particle, double timeElapsed, double epsilon) const;

  Eigen::Vector3d normal;
  double offset, velocity;
};
