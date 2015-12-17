#include "planes.hpp"
#include <iostream>




//outside "semantics" are based on the initial configuration
bool Plane::outside(Particle& particle) const {
   if(particle.restPosition.dot(normal) > offset){
      return true;
   }

   return false;
}

void Plane::doBounce(Particle& p, Eigen::Vector3d n, double w, double epsilon) const {
   p.position += (epsilon + w - p.position.dot(n))*n;

   //zero velocity in the n direction
   //p.velocity -= p.velocity.dot(n)*n;

   //invert velocity in the n direction
   p.velocity -= 2.0*p.velocity.dot(n)*n;

   //friction
   p.velocity *= 0.4;
}


bool Plane::bounceParticle(Particle& particle, double timeElapsed, double epsilon) const {
   bool bounced = false;

   //implementing a double sided bounce first
   if(outside(particle)){
      //check if particle is now INSIDE and then bounce
      if(particle.position.dot(normal) < offset){
         bounced = true;
         //front side bounce
         doBounce(particle, normal, offset, epsilon);
      }
   } else {
      //check if particle is now OUTSIDE and then bounce
      if(particle.position.dot(normal) >= offset){
         bounced = true;
         //back side bounce
         doBounce(particle, -normal, -offset, epsilon);
      }
   }

   return bounced;
}


bool DynamicPlane::bounceParticle(Particle& particle, double timeElapsed, double epsilon) const {
   double currOffset = offset;
   Eigen::Vector3d currNormal = normal;
   bool bounced = false;

   //Only inside particles can bounce
   if(!outside(particle)){
      //check if particle is now OUTSIDE and then bounce
      if(particle.position.dot(currNormal) >= currOffset){
         bounced = true;
         //back side bounce
         doBounce(particle, -currNormal, -currOffset, epsilon);
      }
   }

   return bounced;
}


bool MovingPlane::bounceParticle(Particle& particle, double timeElapsed, double epsilon) const {
   double currOffset = offset + timeElapsed*velocity;
   Eigen::Vector3d currNormal = normal;
   bool bounced = false;

   //Only inside particles can bounce
   if(!outside(particle)){
      //check if particle is now OUTSIDE and then bounce
      if(particle.position.dot(currNormal) >= currOffset){
         bounced = true;
         //back side bounce
         doBounce(particle, -currNormal, -currOffset, epsilon);
      }
   }

   return bounced;
}


void MovingPlane::moveParticle(Particle& particle, double timeElapsed) const {
   if(outside(particle)){
      particle.velocity = velocity*normal;
   }
}


bool TwistingPlane::bounceParticle(Particle& particle, double timeElapsed, double epsilon) const {
   if (timeElapsed > lifetime) return false;

   return DynamicPlane::bounceParticle(particle, timeElapsed, epsilon);
}


void TwistingPlane::twistParticle(Particle& particle, double timeElapsed) const {
   if (timeElapsed > lifetime) return;
   //std::cout<<timeElapsed<<" "<<lifetime<<std::endl;
   if(outside(particle)){
      //will rotate around support point
      const double sos = normal.dot(normal);
      const Eigen::Vector3d supportPoint{normal.x()*offset/sos,
         normal.y()*offset/sos,
         normal.z()*offset/sos};

      //project current position into the plane, compute dir vector from
      //supportPoint to particle position
      Eigen::Vector3d posInPlane = particle.position - supportPoint; 
      posInPlane = posInPlane - particle.position.dot(normal) * normal;

      //compute tangent space, to get velocity components
      Eigen::Vector3d tangent1, tangent2;

      tangent1 = normal.cross(Eigen::Vector3d{1,0,0});
      if(tangent1.isZero(1e-3)){
         tangent1 = normal.cross(Eigen::Vector3d{0,0,1});
         if(tangent1.isZero(1e-3)){
            tangent1 = normal.cross(Eigen::Vector3d{0,1,0});
         }
      }
      tangent1.normalize();

      tangent2 = normal.cross(tangent1);
      tangent2.normalize(); //probably not necessary

      //dot posInPlane with current tangent space to compute rotation velocity
      double t2_rot = angularVelocity * posInPlane.dot(tangent1);
      double t1_rot = angularVelocity * posInPlane.dot(tangent2);

      particle.velocity = -t1_rot * tangent1 + t2_rot * tangent2;
   }
}



bool TiltingPlane::bounceParticle(Particle& particle, double timeElapsed, double epsilon) const {
   bool bounced = false;
   if (timeElapsed > lifetime) return bounced;

   double currOffset = offset;
   Eigen::AngleAxisd t(timeElapsed*angularVelocity,tilt);
   Eigen::Vector3d currNormal = t * normal;

   //Only inside particles can bounce
   if(!outside(particle)){
      //check if particle is now OUTSIDE and then bounce
      if(particle.position.dot(currNormal) >= currOffset){
         bounced = true;
         //back side bounce
         doBounce(particle, -currNormal, -currOffset, epsilon);
      }
   }

   return bounced;
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

