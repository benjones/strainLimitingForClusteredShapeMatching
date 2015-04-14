#include "projectile.hpp"
#include "particle.h"

void Projectile::bounceParticle(Particle& particle, double timeElapsed) const {
  
  Eigen::Vector3d currentPosition = start + timeElapsed*velocity;
  Eigen::Vector3d fromCenter = particle.position - currentPosition;
  if(fromCenter.norm() < radius){
	particle.position = currentPosition + radius*fromCenter.normalized();
	particle.velocity += momentumScale*velocity;
  }

}
