#pragma once

#include "particle.h"
#include "movingPlane.hpp"
#include "twistingPlane.hpp"
#include "tiltingPlane.hpp"
#include "projectile.hpp"
#include "cylinderObstacle.hpp"
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

  void assertFinite() const {
	for(auto& p : particles){
	  assert(p.position.allFinite());
	}
  }


  using FractureInfo = std::tuple<size_t, double,Eigen::Vector3d>;
  void doFracture(std::vector<FractureInfo> potentialSplits);


  inline void restart(){ 
	Eigen::Vector3d oldCameraPosition = cameraPosition;
	Eigen::Vector3d oldCameraLookAt = cameraLookAt;
	Eigen::Vector3d oldCameraUp = cameraUp;

	particles.clear(); 
	planes.clear();
   movingPlanes.clear();
   twistingPlanes.clear();
   tiltingPlanes.clear();
	clusterCenters.clear();
	projectiles.clear();
	cylinders.clear();

	loadFromJson(filename);

	cameraPosition = oldCameraPosition;
	cameraLookAt = oldCameraLookAt;
	cameraUp = oldCameraUp;
	elapsedTime = 0;
  }

  void dumpParticlePositions(const std::string& filename) const;
  void dumpColors(const std::string& filename) const;

  void draw(SDL_Window* window) const ;
  void drawPretty(SDL_Window* window) const ;
  void drawSingleCluster(SDL_Window* window, int frame) const;
  void drawPlanes() const;
  void drawPlanesPretty() const;
  void drawPlane(const Eigen::Vector3d& normal, double offset) const;
  void drawTPlane(const Eigen::Vector3d& normal, double offset, double roffset, double width) const;
  void drawTiltPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& tilt, double offset, double roffset, double width) const;
  void zoom(int amount);
  void pan(Eigen::Vector2i oldposition,
		   Eigen::Vector2i newPosition);

  void move(bool forward);

  void bounceOutOfPlanes();

  void selfCollisions();

  void solveConstraints();
  void updateNeighbors(size_t partIndex);

  void makeClusters();

  void countClusters();

  //pass a container of clusters to update
  template <typename Container>
  void updateClusterProperties(const Container& clusterIndices);

  std::vector<Particle> particles;
  std::vector<Cluster> clusters;

  void printCOM() const;
  
  Eigen::Vector3d computeNeighborhoodCOM(const Cluster& c) const;
  Eigen::Vector3d computeClusterVelocity(const Cluster& c) const;
  
  template <typename cont>
  double sumMass(const cont& indices) const{
	return std::accumulate(indices.begin(), indices.end(),
		0.0,
		[this](double acc, typename cont::value_type index){
		  return acc + particles[index].mass;
		});
  }
  
  template <typename cont> 
  //weighted by # clusters stuff belongs to
  double sumWeightedMass(const cont& indices) const {
	return std::accumulate(indices.begin(), indices.end(), 0.0,
		[this](double acc, typename cont::value_type index){
		  return acc + particles[index].mass/particles[index].numClusters;
		});
  }

  //weighted by # clusters stuff belongs to
  double sumWeightedMass(const std::vector<int> &indices, const std::vector<double> &weights) const {
	double mass = 0.0;
	std::vector<int>::const_iterator ii = indices.begin();
	std::vector<double>::const_iterator di = weights.begin();
	for (; ii!=indices.end(); ii++, di++) mass += (*di)*particles[*ii].mass;
	return mass;
  }

  Eigen::Vector3d sumWeightedRestCOM(const std::vector<int> &indices, const std::vector<double> &weights) const {
	double mass = 0.0;
	Eigen::Vector3d com = Eigen::Vector3d::Zero();
	std::vector<int>::const_iterator ii = indices.begin();
	std::vector<double>::const_iterator di = weights.begin();
	for (; ii!=indices.end(); ii++, di++) {
	  mass += (*di)*particles[*ii].mass;
	  com += (*di)*particles[*ii].restPosition;
	}
	return (com / mass);
  }

  Eigen::Vector3d sumWeightedWorldCOM(const std::vector<int> &indices, const std::vector<double> &weights) const {
	double mass = 0.0;
	Eigen::Vector3d com = Eigen::Vector3d::Zero();
	std::vector<int>::const_iterator ii = indices.begin();
	std::vector<double>::const_iterator di = weights.begin();
	for (; ii!=indices.end(); ii++, di++) {
	  mass += (*di)*particles[*ii].mass;
	  com += (*di)*particles[*ii].position;
	}
	return (com / mass);
  }

  template <typename cont>
  Eigen::Vector3d sumRestCOM(const cont& indices, double totalMass) const{
	return std::accumulate(indices.begin(), indices.end(),
		Eigen::Vector3d::Zero().eval(),
		[this](const Eigen::Vector3d& acc, typename cont::value_type index){
		  return acc + particles[index].mass*
			particles[index].restPosition;
		})/totalMass;
  }

  template <typename cont>
  Eigen::Vector3d sumWeightedRestCOM(const cont& indices, double totalMass) const{
	return 
	  std::accumulate(indices.begin(), indices.end(),
		  Eigen::Vector3d::Zero().eval(),
		  [this](const Eigen::Vector3d& acc, typename cont::value_type index){
			return acc + (particles[index].mass/particles[index].numClusters)*
			  particles[index].restPosition;
		  }) /
	  totalMass;
  }



  template <typename cont>
  Eigen::Vector3d sumWorldCOM(const cont& indices, double totalMass) const{
	return std::accumulate(indices.begin(), indices.end(),
		Eigen::Vector3d::Zero().eval(),
		[this](const Eigen::Vector3d& acc, typename cont::value_type index){
		  return acc + particles[index].mass*
			particles[index].position;
		})/totalMass;
  }


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
  std::vector<MovingPlane> movingPlanes;
  std::vector<TwistingPlane> twistingPlanes;
  std::vector<TiltingPlane> tiltingPlanes;
  std::vector<Projectile> projectiles;
  std::vector<CylinderObstacle> cylinders;
  std::vector<size_t> clusterCenters;

  double dt, elapsedTime;
  double neighborRadius;
  int nClusters;

  int numConstraintIters;
  double omega, gamma, alpha, springDamping;
  double yield, nu, hardening; // plasticity parameters
  double toughness;
  int maxNumClusters;

  bool drawClusters = true;
  bool drawColoredParticles = false;
  bool colorByToughness = false;
  bool dragWithPlanes = true;
  bool paused = false;
  int which_cluster = -1;

  benlib::Profiler prof;
};
