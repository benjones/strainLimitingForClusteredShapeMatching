#pragma once
#include <Eigen/Dense>
#include "particle.h"

//static plane class -- bounces on both sides
class Plane {
   public:
      Plane(Eigen::Vector3d _normal, double _offset) :
         normal(_normal), offset(_offset) {}

      virtual bool outside(Particle& particle) const;
      void doBounce(Particle& p, Eigen::Vector3d n, double w, double epsilon) const;
      virtual bool bounceParticle(Particle& particle, double timeElapsed, double epsilon) const;

      Eigen::Vector3d normal;
      double offset;
};

//dynamic planes affect the dynamics of particles on their outside
//bounce only on inside
class DynamicPlane : public Plane {
   public:
      DynamicPlane(Eigen::Vector3d _normal, double _offset) :
         Plane(_normal, _offset) {}

      virtual bool bounceParticle(Particle& particle, double timeElapsed, double epsilon) const;
};


class TwistingPlane : public DynamicPlane {
   public:
      TwistingPlane(Eigen::Vector3d _normal, double _offset, double _angularVelocity, double _width, double _lifetime)
         : DynamicPlane(_normal, _offset), angularVelocity(_angularVelocity), width(_width), lifetime(_lifetime)
      {}

      virtual bool bounceParticle(Particle& particle, double timeElapsed, double epsilon) const;
      void twistParticle(Particle& particle, double timeElapsed) const;

      double angularVelocity;
      double width;
      double lifetime;
};


class TiltingPlane : public DynamicPlane {
   public:

      TiltingPlane(Eigen::Vector3d _normal, Eigen::Vector3d _tilt, double _offset, double _angularVelocity, double _width, double _lifetime)
         : DynamicPlane(_normal, _offset), tilt(_tilt.normalized()), angularVelocity(_angularVelocity), width(_width), lifetime(_lifetime)
      {}

      virtual bool bounceParticle(Particle& particle, double timeElapsed, double epsilon) const;
      void tiltParticle(Particle& particle, double timeElapsed) const;

      Eigen::Vector3d tilt;
      double angularVelocity;
      double width;
      double lifetime;
};

class MovingPlane : public DynamicPlane {
   public:
      MovingPlane(Eigen::Vector3d _normal, double _offset, double _velocity)
         : DynamicPlane(_normal, _offset), velocity(_velocity)
      {}

      virtual bool bounceParticle(Particle& particle, double timeElapsed, double epsilon) const;
      void moveParticle(Particle& particle, double timeElapsed) const;

      double velocity;
};
