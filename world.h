#pragma once

#include <iostream>
#include "particle.h"
#include "movingPlane.hpp"
#include "twistingPlane.hpp"
#include "tiltingPlane.hpp"
#include "projectile.hpp"
#include "cylinderObstacle.hpp"

#include "preallocVector.hpp"
#include "profiler.hpp"
#include <SDL.h>

#include "clustering.h"

class World{
public:
  std::string filename; //for reloading the file

  void loadFromJson(const std::string& _filename);
  

  struct ParticleSet{ 
	std::vector<Eigen::Vector3d> positions;
	Eigen::Vector3d bbMin, bbMax;
  };

  ParticleSet readParticleFile(const std::string& _filename);
  void saveParticleFile(const std::string& _filename) const;
  
  void timestep();

  void assertFinite() const {
	for(auto& p : particles){
	  assert(p.position.allFinite());
	}
  }


  struct FractureInfo{
	size_t clusterIndex;
	double effectiveToughness;
	Eigen::Vector3d fractureNormal;
	 
  };
  void doFracture(std::vector<FractureInfo> potentialSplits);
  void splitOutliers();
  void cullSmallClusters();
  void removeLonelyParticles();

  inline void restart(){ 
	/*	Eigen::Vector3d oldCameraPosition = cameraPosition;
	Eigen::Vector3d oldCameraLookAt = cameraLookAt;
	Eigen::Vector3d oldCameraUp = cameraUp;
	*/
	particles.clear(); 
	clusters.clear();
	planes.clear();
	movingPlanes.clear();
	twistingPlanes.clear();
	tiltingPlanes.clear();
	projectiles.clear();
	cylinders.clear();
	clusterCenters.clear();

	loadFromJson(filename);
	initializeNeighbors();
	/*
	cameraPosition = oldCameraPosition;
	cameraLookAt = oldCameraLookAt;
	cameraUp = oldCameraUp;
	*/
	prof.dump<std::chrono::duration<double>>(std::cout);
	prof.reset();
	elapsedTime = 0;
  }

  void dumpParticlePositions(const std::string& filename) const;
  void dumpClippedSpheres(const std::string& filename) const;
  void dumpColors(const std::string& filename) const;


  //refactored out into vis, camera/settings stuff
  /*  void draw(SDL_Window* window) const ;
  void drawPretty(SDL_Window* window) const ;
  void drawSingleCluster(SDL_Window* window, int frame) const;
  void drawPlanes() const;
  void drawPlanesPretty() const;
  void drawTPlane(const Eigen::Vector3d& normal, double offset, double roffset, double width) const;
  void drawTiltPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& tilt, double offset, double roffset, double width) const;clus
  void zoom(int amount);
  void pan(Eigen::Vector2i oldposition,
		   Eigen::Vector2i newPosition);

  void move(bool forward);
  */
  void bounceOutOfPlanes();

 // std::vector<std::vector<int> > clusterCollisionMap;
 //std::unordered_set<std::pair<int,int>> clusterCollisionMap;
  //void buildClusterMaps();
  void initializeNeighbors();
  void selfCollisions();

  void solveConstraints();
  void updateNeighbors(size_t partIndex);


  void setupPlaneConstraints();

  void countClusters();
  void mergeClusters(const std::vector<Particle>& newParticles,
	  const std::vector<Cluster>& newClusters);

  //pass a container of clusters to update
  template <typename Container>
  void updateClusterProperties(const Container& clusterIndices);

  std::vector<Particle> particles;
  std::vector<Cluster> clusters;

  void printCOM() const;
  
  Eigen::Vector3d computeClusterVelocity(const Cluster& c) const;
  
  template <typename cont>
  double sumMass(const cont& indices) const{
	return std::accumulate(indices.begin(), indices.end(),
		0.0,
		[this](double acc, typename cont::value_type index){
		  return acc + particles[index].mass;
		});
  }

  double sumWeightedMass(const std::vector<Cluster::Member > &members) const {
	double mass = 0.0;
	for (auto &member : members) {
	  auto &p = particles[member.index];
	  mass += (member.weight / p.totalweight) * p.mass;
	}
	return mass;
  }

  Eigen::Vector3d sumWeightedRestCOM(const std::vector<Cluster::Member> &members) const {
	double mass = 0.0;
	Eigen::Vector3d com = Eigen::Vector3d::Zero();
	for (auto &member : members) {
	  auto &p = particles[member.index];
	  double w = (member.weight / p.totalweight) * p.mass;
	  mass += w;
	  com += w * p.restPosition;
	}
	return (com / mass);
  }

  Eigen::Vector3d sumWeightedWorldCOM(const std::vector<Cluster::Member>& members) const {
	double mass = 0.0;
	Eigen::Vector3d com = Eigen::Vector3d::Zero();
	for (auto &member : members) {
	  auto &p = particles[member.index];
	  double w = (member.weight / p.totalweight) * p.mass;
	  mass += w;
	  com += w * p.position;
	}
	return (com / mass);
  }

  //compute the APQ matrix (see meuller 2005)
  Eigen::Matrix3d computeApq(const Cluster& c) const;
  void updateTransforms(Cluster& c) const;
  
  inline Eigen::Vector3d getMomentum(){
	return std::accumulate(particles.begin(), particles.end(),
						   Eigen::Vector3d{0.0,0.0,0.0},
						   [](Eigen::Vector3d acc, const Particle& p){
							 return acc + p.velocity;
						   });
  }


  void strainLimitingIteration();



  Eigen::Vector3d gravity;

  std::vector<Eigen::Vector4d> planes;
  std::vector<MovingPlane> movingPlanes;
  std::vector<TwistingPlane> twistingPlanes;
  std::vector<TiltingPlane> tiltingPlanes;
  std::vector<Projectile> projectiles;
  std::vector<CylinderObstacle> cylinders;
  std::vector<size_t> clusterCenters;
  bool dragWithPlanes = true;



  double dt, elapsedTime;
  double neighborRadius, neighborRadiusMax;

  ClusteringParams clusteringParams;

  bool fractureOn, selfCollisionsOn, delayRepeatedFracture;

  
  int numConstraintIters;
  double omega, gamma, alpha, springDamping;
  double yield, nu, hardening; // plasticity parameters
  double toughness, toughnessBoost, toughnessFalloff;
  double collisionRestitution;
  double collisionGeometryThreshold;
  double outlierThreshold;



  benlib::Profiler prof;
};
