#include "cylinderObstacle.hpp"
#include "particle.h"

void CylinderObstacle::bounceParticle(Particle& particle) const {
  Eigen::Vector3d projectedPoint =
	supportPoint + (particle.position - supportPoint).dot(normal)*normal;
  
  Eigen::Vector3d fromProjected = particle.position - projectedPoint;
  if( fromProjected.norm() < radius){
	//bounce
	particle.position = projectedPoint + radius*fromProjected.normalized();
  }
  
}
