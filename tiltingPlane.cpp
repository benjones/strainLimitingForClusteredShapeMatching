#include "tiltingPlane.hpp"
#include <iostream>


bool TiltingPlane::outside(Particle& particle) const {
   if(particle.restPosition.dot(normal) > offset){
      return true;
   }

   return false;
}


void TiltingPlane::tiltParticle(Particle& particle, double timeElapsed) const {
   if (timeElapsed > lifetime) return;
   //std::cout<<timeElapsed<<" "<<lifetime<<std::endl;
   if(outside(particle)){
      //will rotate around support point
      const double sos = normal.dot(normal);
      const Eigen::Vector3d supportPoint{normal.x()*offset/sos,
         normal.y()*offset/sos,
         normal.z()*offset/sos};

      //compute tangent space, to get velocity components
      Eigen::Vector3d tangent;

      tangent = normal.cross(tilt);

      //project current position into the plane defined by tilt, compute dir vector from
      //supportPoint to particle position
      Eigen::Vector3d posInPlane = particle.position - supportPoint; 
      posInPlane = posInPlane - particle.position.dot(tilt) * tilt;

      //dot posInPlane with current tangent space to compute rotation velocity
      double t2_rot = angularVelocity * posInPlane.dot(normal);
      double t1_rot = angularVelocity * posInPlane.dot(tangent);
      
      particle.velocity = t1_rot * normal - t2_rot * tangent;
   }
}

