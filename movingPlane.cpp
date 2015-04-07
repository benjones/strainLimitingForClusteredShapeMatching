#include "movingPlane.hpp"
#include "particle.h"


void MovingPlane::bounceParticle(Particle& particle, double timeElapsed) const {

  double currentOffset = offset + timeElapsed*velocity;
  
  if(particle.restPosition.dot(normal) > offset){
	
	Eigen::Vector3d tangential = 
	  particle.restPosition - (particle.restPosition.dot(normal))*normal;

	particle.position = tangential + currentOffset*normal;
	
  }

}
