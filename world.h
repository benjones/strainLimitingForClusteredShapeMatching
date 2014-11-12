#pragma once

#include "particle.h"
#include "accelerationGrid.h"
#include "preallocVector.hpp"
#include "profiler.hpp"
#include <SDL.h>

class World{
public:
  std::string filename; //for reloading the file

  void loadFromJson(const std::string& _filename);
  
  void readParticleFile(const std::string& _filename);
  void saveParticleFile(const std::string& _filename) const;
  
  void timestep();
  inline void restart(){ 
	Eigen::Vector3d oldCameraPosition = cameraPosition;
	Eigen::Vector3d oldCameraLookAt = cameraLookAt;
	Eigen::Vector3d oldCameraUp = cameraUp;

	particles.clear(); 
	planes.clear();
	clusterCenters.clear();
	
	loadFromJson(filename);

	cameraPosition = oldCameraPosition;
	cameraLookAt = oldCameraLookAt;
	cameraUp = oldCameraUp;
  }

  void draw(SDL_Window* window) const ;
  void zoom(int amount);
  void pan(Eigen::Vector2i oldposition,
		   Eigen::Vector2i newPosition);

  void move(bool forward);

  void bounceOutOfPlanes();

  void solveConstraints();
  void updateNeighbors(size_t partIndex);

  void makeClusters();

  std::vector<Particle> particles;
  std::vector<Cluster> clusters;

  void printCOM() const;
  
  Eigen::Vector3d computeNeighborhoodCOM(const Cluster& c) const;
  Eigen::Vector3d computeClusterVelocity(const Cluster& c) const;
  
  //compute the APQ matrix (see meuller 2005)
  Eigen::Matrix3d computeApq(const Cluster& c, 
							 const Eigen::Matrix3d& init,
							 const Eigen::Vector3d& worldCOM) const;
  
  inline Eigen::Vector3d getMomentum(){
	return std::accumulate(particles.begin(), particles.end(),
						   Eigen::Vector3d{0.0,0.0,0.0},
						   [](Eigen::Vector3d acc, const Particle& p){
							 return acc + p.velocity;
						   });
  }


  void strainLimitingIteration();



  struct PositionGetter{
	Eigen::Vector3d operator()(const Particle& p){
	  return p.position;
	}
  };
  
  struct RestPositionGetter{
	Eigen::Vector3d operator()(const Particle& p){
	  return p.position;
	}
  };
  AccelerationGrid<Particle, PositionGetter> positionGrid;
  AccelerationGrid<Particle, RestPositionGetter> restPositionGrid;
  

  Eigen::Vector3d cameraPosition, cameraLookAt, cameraUp, gravity;

  std::vector<Eigen::Vector4d> planes;

  std::vector<size_t> clusterCenters;

  double dt;
  double neighborRadius;
  int nClusters;

  int numConstraintIters;
  double omega, gamma, alpha, springDamping;

  int maxNumClusters;

  benlib::Profiler prof;
};
