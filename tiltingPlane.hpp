#pragma once
#include <Eigen/Dense>
#include "particle.h"

class TiltingPlane{
public:

  TiltingPlane(Eigen::Vector3d _normal, Eigen::Vector3d _tilt, double _offset, double _angularVelocity, double _width, double _lifetime)
	:normal(_normal.normalized()), tilt(_tilt.normalized()), offset(_offset), angularVelocity(_angularVelocity), width(_width), lifetime(_lifetime)
  {}

  bool outside(Particle& particle) const;
  void tiltParticle(Particle& particle, double timeElapsed) const;

  Eigen::Vector3d normal;
  Eigen::Vector3d tilt;
  double offset;
  double angularVelocity;
  double width;
  double lifetime;
};
