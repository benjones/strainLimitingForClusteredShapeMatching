#pragma once
#include <Eigen/Dense>
#include <vector>
#include <unordered_set>

#include "color_spaces.h"

namespace Ogre{
  class SceneManager;
  class SceneNode;
  class Entity;
  class Material;
}

class Particle;
inline void noop(Particle&){}

class Particle {
public:
  Eigen::Vector3d position, 
	velocity, 
	restPosition,  // DO NOT USE for anything other than moving planes
	oldPosition, 
	goalPosition, 
	goalVelocity,
	embeddedPosition;
  double mass; 
  double totalweight; // sum of weights of this particle in all clusters
  double tmpd;
  size_t numClusters; // only used to determine singelton particles
  bool outsideSomeMovingPlane;
  std::vector<int> clusters; 
  Particle() {};
  Particle(const Particle &p)  = default; //c++11 magic :)
  short flags;
  static const short SPLIT = 1;
  static const short JUST_SPLIT = 2;
  int id;
  Ogre::SceneNode* sceneNode = nullptr;
  Ogre::Entity* entity = nullptr;
  //Ogre::Material* material = nullptr;
  //cleanup stuff without knowing anything about ogre
  std::function<void(Particle&)> cleanup = noop;
  RGBColor color;
  
};

class CollisionGeometry {
 private:
 public:
  Eigen::Vector3d c;
  double r;
  CollisionGeometry(){};
  CollisionGeometry & operator= (const CollisionGeometry &that);
  CollisionGeometry(const CollisionGeometry &that) = default;
  bool project(Eigen::Vector3d &x);
  void init (const Eigen::Vector3d &c, double r) { this->c = c; this->r = r;};
  void addPlane(const Eigen::Vector3d &n, double offset) { planes.emplace_back(n, offset); };
  
  std::vector<std::pair<Eigen::Vector3d, double> > planes;
};

class Cluster {
 public:
  Eigen::Vector3d restCom; // DO NOT USE except for initial clustering
  Eigen::Vector3d worldCom;
  Eigen::Matrix3d aInv, Fp, FpNew;
  //contains the particle index and its weight
  struct Member{
  	int index; double weight;  Eigen::Vector3d pos;
  };
  std::vector<Member > members;
  std::unordered_set<int> neighbors;  // Note: this refers to neighboring CLUSTERS 
  std::unordered_set<int> oldNeighbors;		// Note: this refers to neighboring CLUSTERS
  double mass, width, renderWidth;
  double toughness;
  double cstrain; // cumulative strain (for work hardening)
  Cluster() {Fp.setIdentity(); cstrain = 0.0;}
  Cluster(const Cluster &c) = default;
  void updatePlasticity(const Eigen::Vector3d &sigma, const Eigen::Matrix3d &U, const Eigen::Matrix3d &V, 
	  double yield, double nu, double hardening);
  
  
  CollisionGeometry cg;
  Eigen::Matrix3d worldToRestTransform, restToWorldTransform;
  Eigen::Matrix4d getVisTransform() const;
  bool justFractured = false;
  double timeSinceLastFracture;
  bool markedForRemoval = false;
  bool newCluster = false;
  int fadeSteps;
  std::vector<double> oweights;
  int embedId;
};

