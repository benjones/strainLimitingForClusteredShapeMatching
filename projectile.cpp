#include "projectile.hpp"
#include "particle.h"

void Projectile::bounceParticle(Particle& particle, double timeElapsed) const {
  
  Eigen::Vector3d currentPosition = start + timeElapsed*velocity;
  Eigen::Vector3d fromCenter = particle.position - currentPosition;
  if(fromCenter.norm() < radius){
	particle.position = currentPosition + 1.1*radius*fromCenter.normalized();
	particle.velocity += 10.0*velocity;
  }

}
