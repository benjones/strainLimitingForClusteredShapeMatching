#pragma once
#include <Eigen/Dense>
#include "particle.h"

class TwistingPlane{
public:

  TwistingPlane(Eigen::Vector3d _normal, double _offset, double _angularVelocity, double _width)
	:normal(_normal.normalized()), offset(_offset), angularVelocity(_angularVelocity), width(_width)
  {}

  bool outside(Particle& particle) const;
  void twistParticle(Particle& particle, double timeElapsed) const;
  void backsideReflectBounceParticle(Particle& particle, double timeElapsed, double epsilon) const;

  Eigen::Vector3d normal;
  double offset;
  double angularVelocity;
  double width;
};
