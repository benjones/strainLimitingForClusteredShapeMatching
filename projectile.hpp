#pragma once

#include <Eigen/Dense>

class Particle;

class Projectile {

public:

  Eigen::Vector3d start, velocity;
  double radius, momentumScale;  
  Projectile(const Eigen::Vector3d& _start, const Eigen::Vector3d& _velocity, 
	  double _radius, double _momentumScale)
	:start{_start}, velocity{_velocity}, radius{_radius}, momentumScale{_momentumScale}
  {}

  void bounceParticle(Particle& particle, double timeElapsed) const;

};
