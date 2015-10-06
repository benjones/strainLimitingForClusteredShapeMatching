#pragma once

#include "world.h"
#include "accelerationGrid.h"
#include "particle.h"

struct RestPositionGetter{
  Eigen::Vector3d operator()(const Particle& p) const{
	return p.position;
  }
};

struct ClusteringParams{

};

std::vector<Cluster> makeClusters(std::vector<Particle>& particles, const ClusteringParams& params);

//Will mody particleSet in place to set weights, and cluster
std::vector<Cluster> makeRandomClusters(std::vector<Particle>& particlest, double neighborRadius);
