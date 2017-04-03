#include "constraintPlane.hpp"
#include "particle.h"
#include "enumerate.hpp"

ConstraintPlane::ConstraintPlane(double xValue, const std::vector<Particle>& particles){

  for(auto&& pr : benlib::enumerate(particles)){
	if(pr.second.position.x() < xValue){
	  constraints.push_back({pr.first, pr.second.position});
	}
  }
}


void ConstraintPlane::constrainParticles( std::vector<Particle>& particles) const{
  for(auto& pr : constraints){
	particles[pr.first].position = pr.second;
	particles[pr.first].velocity.setZero();
  }
}
