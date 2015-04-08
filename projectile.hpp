#pragma once

#include <Eigen/Dense>

class Particle;

class Projectile {

public:

  Eigen::Vector3d start, velocity;
  double radius;  
  Projectile(const Eigen::Vector3d& _start, const Eigen::Vector3d& _velocity, double _radius)
	:start{_start}, velocity{_velocity}, radius{_radius}
  {}

  void bounceParticle(Particle& particle, double timeElapsed) const;

};
