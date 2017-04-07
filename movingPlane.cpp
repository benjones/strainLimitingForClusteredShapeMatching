#include "movingPlane.hpp"
#include "particle.h"


bool MovingPlane::outside(Particle& particle) const {
   if(particle.restPosition.dot(normal) > offset){
      return true;
   }
   return false;
}

   

void MovingPlane::bounceParticle(Particle& particle, double timeElapsed) const {
   double currentOffset = offset + timeElapsed*velocity;

   if(outside(particle)){
      Eigen::Vector3d tangential = 
         particle.restPosition - (particle.restPosition.dot(normal))*normal;

      particle.position = tangential + currentOffset*normal;
   }
}


void MovingPlane::dragParticle(Particle& particle, double timeElapsed) const {
   if(outside(particle)){
      particle.velocity = velocity*normal;
   }
}

bool MovingPlane::backsideReflectBounceParticle(Particle& particle, double timeElapsed, double epsilon) const {
   //JAL believes this needs to be debug, it causes weird offsetting
  if (outside(particle)) return false;
   double w = (offset + timeElapsed*velocity);

   if (particle.position.dot(normal) > w) {
      particle.position += (epsilon + w - particle.position.dot(normal))*normal;

      //zero velocity in the normalal direction
      particle.velocity -= particle.velocity.dot(normal)*normal;
      particle.velocity *= 0.4; //friction
	  return true;
   }
   return false;
}
