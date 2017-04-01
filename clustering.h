#pragma once

#include "particle.h"

struct RestPositionGetter{
  Eigen::Vector3d operator()(const Particle& p) const{
	return p.restPosition;
  }
};

struct EmbeddedPositionGetter{
  Eigen::Vector3d operator()(const Particle& p) const{
	return p.embeddedPosition;
  }
};

struct ClusteringParams{
  double neighborRadius;
  int nClusters;
  double neighborRadiusMax;
  int nClustersMax;
  int clusterItersMax;
  int clusteringAlgorithm; // 0 = default (fuzzy c-means with weights); 1 = k-means; 2 = random
  double clusterOverlap;
  int clusterKernel; // 0 = 1 / (r^2 + eps), 1 = constant, 2 = poly6,
                     //3 = constant/poly6 blend (with kernelweight), 4 = fuzzy c-means

  double kernelWeight;// only for the blended kernel, and fuzzy c-means
  double blackhole;
  double poly6norm, sqrNeighborRadius;

  double kernel(const Eigen::Vector3d &x) const;

};

//params will change if clustering fails
std::vector<Cluster> iterateMakeClusters(std::vector<Particle>& particles,  ClusteringParams& params);

std::vector<Cluster> makeClusters(std::vector<Particle>& particles, const ClusteringParams& params);

//Will mody particleSet in place to set weights, and cluster
std::vector<Cluster> makeRandomClusters(std::vector<Particle>& particles, double neighborRadius);
