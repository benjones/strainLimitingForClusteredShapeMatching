#include "clustering.h"
#include "range.hpp"
#include "utils.h"

std::vector<Cluster> makeRandomClusters(std::vector<Particle>& particle, double neighborRadius) {
  std::random_device rd;
  std::mt19937 g(std::mt19937::default_seed);

  std::vector<Cluster> clusters;
  
  auto r = benlib::range(particles.size());
  auto lonelyParticles = std::vector<size_t>(r.begin(), r.end());

  AccelerationGrid<Particle, RestPositionGetter> restPositionGrid;
  restPositionGrid.numBuckets = 8; //TODO tune me, 512 buckets now... seems reasonable?
  restPositionGrid.updateGrid(particles);

  
  for(auto& p : particles){p.numClusters = 0; p.totalweight = 0.0;}

  
  while(!lonelyParticles.empty()) {
	int r = g() % lonelyParticles.size();
	//auto currentParticle = lonelyParticles.back();
	//lonelyParticles.pop_back();
	auto currentParticle = lonelyParticles[r];
	lonelyParticles.erase(lonelyParticles.begin() + r);
	Cluster c;
	c.restCom = particles[currentParticle].restPosition;
	std::vector<int> neighbors = restPositionGrid.getNearestNeighbors(particles, c.restCom, neighborRadius);
	c.neighbors.resize(neighbors.size());
	double w = 1.0;
	for(unsigned int i=0; i<neighbors.size(); i++) {
	  auto n = neighbors[i];
	  Particle &p = particles[n];
	  p.numClusters++; 
	  p.totalweight += w; 
	  p.clusters.push_back(clusters.size());
	  c.neighbors[i] = std::pair<int, double>(n, w);
	}
	clusters.push_back(c);

	//std::cout<<lonelyParticles.size()<<std::endl;
	auto it = std::remove_if(lonelyParticles.begin(), lonelyParticles.end(),
		[&neighbors](size_t n){
		  //neighbors aren't lonely anymore
		  return utils::containsValue(neighbors, n);
		});
	lonelyParticles.erase(it, lonelyParticles.end());
	//std::cout<<lonelyParticles.size()<<std::endl;
  } 
  for (auto& c : clusters) c.cg.init(c.restCom, neighborRadius);
  for (auto& c : clusters) {
	c.mass = 0.0;
	c.restCom = Eigen::Vector3d::Zero();
	for (auto &n : c.neighbors) {
	  auto &p = particles[n.first];
	  double w = n.second / p.totalweight;
	  c.restCom += w * p.mass * p.position;
	  c.mass += w * p.mass;
	}
	c.restCom /= c.mass;
  }
  return true;
}

std::vector<Cluster> World::makeClusters(World::ParticleSet& particleSet,
	const ClusteringParams& params){
  
  std::cout<<"using clusteringAlgorithm "<<clusteringAlgorithm<<std::endl;
  std::vector<Cluster> clusters;

  bool converged;
  int iters = 0;
  double convergenceThreshold = 1e-6 * sqrNeighborRadius; // 0.1% motion allowed
  std::random_device rd;
  std::mt19937 g(std::mt19937::default_seed);
  std::vector<bool> picked(particles.size(), false);

  // initialize cluster centers to random paricles
  if (clusteringAlgorithm == 2) {
	makeRandomClusters();
	converged = true;
  }

  
  if (clusteringAlgorithm < 2) {
	while (clusters.size() < nClusters) {
	  Cluster c;
	  int r = g() % particles.size();
	  if (picked[r]) continue;
	  picked[r] = true;
	  c.restCom = particles[r].restPosition;
	  clusters.push_back(c);
	}

	converged = false;
	while (!converged) {
	  iters++;
	  converged = true;
	  for (auto& c : clusters) c.neighbors.clear();
	  
	  for (auto i=0; i<particles.size(); i++) {
		auto& p = particles[i];
		int bestCluster = 0;
		double bestNorm = (clusters[0].restCom - p.restPosition).squaredNorm();
		for (auto j = 0; j<clusters.size(); j++) {
		  double newNorm = (clusters[j].restCom - p.restPosition).squaredNorm();
		  if (newNorm < bestNorm) {
			bestCluster = j;
			bestNorm = newNorm;
		  }
		}
		clusters[bestCluster].neighbors.push_back(std::pair<int,double>(i,1.0));
		if (p.numClusters != bestCluster) converged = false;
		p.numClusters = bestCluster;
		p.totalweight = 1.0;
	  }
	  
	  for (auto& c : clusters) {
		c.mass = 0.0;
		c.restCom = Eigen::Vector3d::Zero();
		for (auto &n : c.neighbors) {
		  auto &p = particles[n.first];
		  double w = n.second / p.totalweight;
		  c.restCom += w * p.mass * p.position;
		  c.mass += w * p.mass;
		}
		c.restCom /= c.mass;
	  }
	}
  }
  
  //double sqrNeighborRadius = neighborRadius*neighborRadius;
  if (clusteringAlgorithm == 1) {
	for(auto& p : particles){p.numClusters = 0; p.totalweight = 0.0;}
	for (auto j=0; j<clusters.size(); j++) {
	  auto &c = clusters[j];
	  std::vector<int> neighbors = restPositionGrid.getNearestNeighbors(particles, c.restCom, neighborRadius);
	  c.neighbors.resize(neighbors.size());
	  for(unsigned int i=0; i<neighbors.size(); i++) {
		auto n = neighbors[i];
		Particle &p = particles[n];
		//double norm = (c.restCom - p.restPosition).squaredNorm();
		//double w = 1e7*cube(sqrNeighborRadius-norm);
		double w = kernel(c.restCom - p.restPosition); //1.0;
		c.neighbors[i] = std::pair<int, double>(n, w);
		p.totalweight += w;
		p.clusters.push_back(j);
		p.numClusters++;
	  }
	}

	for (auto& c : clusters) {
	  c.cg.init(c.restCom, neighborRadius);
	  c.mass = 0.0;
	  c.restCom = Eigen::Vector3d::Zero();
	  for (auto &n : c.neighbors) {
		auto &p = particles[n.first];
		double w = n.second / p.totalweight;
		c.restCom += w * p.mass * p.position;
		c.mass += w * p.mass;
	  }
	  c.restCom /= c.mass;
	}
	for (auto& p : particles) {
	  if (p.numClusters == 0) {
		std::cout<<"Particle has no cluster"<<std::endl;
		return false;
	  }
	}
  }

  // fuzzy c-means loop
  if (clusteringAlgorithm < 1) {
	iters = 0;
	while ((!converged || iters < 5) && iters < clusterItersMax) {
	  converged = true;
	  iters++;
	
	  for (auto i=0; i<particles.size(); i++) {
		auto& p = particles[i];
		p.clusters.clear();
		p.numClusters = 0;
		p.totalweight = 0.0;
	  }
	  
	  for (auto j = 0; j<clusters.size(); j++) {
		Cluster &c = clusters[j];
		int oldNeighbors = c.neighbors.size();
		std::vector<int> neighbors = restPositionGrid.getNearestNeighbors(particles, c.restCom, neighborRadius);
		c.neighbors.resize(neighbors.size());
		for(unsigned int i=0; i<neighbors.size(); i++) {
		  auto n = neighbors[i];
		  Particle &p = particles[n];
		  double w = kernel (c.restCom - p.restPosition);
		  c.neighbors[i] = std::pair<int, double>(n, w);
		  p.totalweight += w;
		  p.clusters.push_back(j);
		  p.numClusters++;
		}
		if (oldNeighbors != c.neighbors.size()) {
		  //std::cout<<oldNeighbors<<" != "<<c.neighbors.size()<<std::endl;
		  converged = false;
		}
	  }
	
	  for (int i = 0; i<particles.size(); i++) {
		auto& p = particles[i];
		if (p.numClusters > 0) continue;
		int bestCluster = 0;
		double bestNorm = (clusters[0].restCom - p.restPosition).squaredNorm();
		for (auto j = 0; j<clusters.size(); j++) {
		  double newNorm = (clusters[j].restCom - p.restPosition).squaredNorm();
		  if (newNorm < bestNorm) {
			bestCluster = j;
			bestNorm = newNorm;
		  }
		}
		clusters[bestCluster].neighbors.push_back(std::pair<int,double>(i,blackhole));
		p.totalweight = 1.0;
		//std::cout<<"particle "<<i<<" not in cluster"<<std::endl;
		converged = false;
	  }

	  if (clusterKernel == 4) {
		for (auto m = 0; m < particles.size(); m++) {
		  auto &p = particles[m];
		  p.totalweight = 0.0;
		  std::vector<double> updatedweights;
		  updatedweights.resize(p.clusters.size());

		  for (auto i=0; i<p.clusters.size(); i++) {
			double sum = 0.0;
			int k = 0;
			for (;k < clusters[p.clusters[i]].neighbors.size() && clusters[p.clusters[i]].neighbors[k].first != m; k++);
			for (auto j=0; j<p.clusters.size(); j++) {
			  int l = 0;
			  for (;l < clusters[p.clusters[j]].neighbors.size() && clusters[p.clusters[j]].neighbors[l].first != m; l++);
			  sum += pow(clusters[p.clusters[i]].neighbors[k].second / clusters[p.clusters[j]].neighbors[l].second, 2.0/kernelWeight-1.0);
			}
			updatedweights[i] = 1.0/sum;
		  }
		  for (auto i=0; i<p.clusters.size(); i++) {
			int k = 0;
			for (;k < clusters[p.clusters[i]].neighbors.size() && clusters[p.clusters[i]].neighbors[k].first != m; k++);
			clusters[p.clusters[i]].neighbors[k].second = updatedweights[i];
		  }
		}
		for (auto &c : clusters) {
		  for (auto &n : c.neighbors) {
			particles[n.first].totalweight += n.second;
		  }
		}
	  }
	  
	  for (auto& c : clusters) {
		c.worldCom = c.restCom;
		c.mass = 0.0;
		c.restCom = Eigen::Vector3d::Zero();
		for (auto &n : c.neighbors) {
		  auto &p = particles[n.first];
		  //double w = n.second / p.totalweight + clusterOverlap / p.numClusters; //* std::max(1.0, p.numClusters * clusterOverlap);
		  double w = n.second / p.totalweight;
		  c.restCom += w * p.mass * p.position;
		  c.mass += w * p.mass;
		}
		c.restCom /= c.mass;
		if ((c.restCom - c.worldCom).squaredNorm() > convergenceThreshold) {
		  //std::cout<<c.restCom(0)<<" "<<c.restCom(1)<<" "<<c.restCom(2)<<" "<<c.worldCom(0)<<" "<<c.worldCom(1)<<" "<<c.worldCom(2)<<" "<<convergenceThreshold<<" "<<(c.restCom - c.worldCom).squaredNorm()<<std::endl;
		  converged = false;
		}
	  }
	} 
	for (auto& c : clusters) c.cg.init(c.restCom, neighborRadius);
	std::cout<<"Fuzzy c-means clustering ran "<<iters<<" iterations."<<std::endl;
  }

  for (auto& p : particles) {
	if (p.numClusters == 0) {
	  std::cout<<"Particle has no cluster"<<std::endl;
	  return false;
	}
  }

  for (auto& c : clusters) {
	Eigen::Matrix3d Aqq;
	Aqq.setZero();
	for (auto &n : c.neighbors) {
	  auto &p = particles[n.first];
	  Eigen::Vector3d qj = p.restPosition - c.restCom;
	  Aqq += qj * qj.transpose();
	}
	Eigen::JacobiSVD<Eigen::Matrix3d> solver(Aqq, Eigen::ComputeFullU | Eigen::ComputeFullV);

	double min, max;
	Eigen::Vector3d n;
	for (unsigned int i=0; i<3; i++) {
	  n = solver.matrixV().col(i);
	  min = std::numeric_limits<double>::max();
	  max = -std::numeric_limits<double>::max();
	  for (auto &neighbor : c.neighbors) {
		auto &p = particles[neighbor.first];
		double d = n.dot(p.restPosition);
		if (d < min) min = d;
		if (d > max) max = d;
	  }
	  double center = n.dot(c.cg.c);
	  //std::cout<<max - min<<" "<<neighborRadius<<std::endl;
	  //if (max - min < 1.5*neighborRadius) {
	  if (max - center < collisionGeometryThreshold*neighborRadius) {
		std::cout<<" adding plane: "<<n(0)<<"x + "<<n(1)<<"y + "<<n(2)<<"z + "<<-max<<std::endl;
		c.cg.addPlane(n, -max);
	  }
	  if (center - min < collisionGeometryThreshold*neighborRadius) {
		std::cout<<" adding plane: "<<-n(0)<<"x + "<<-n(1)<<"y + "<<-n(2)<<"z + "<<min<<std::endl;
		c.cg.addPlane(-n, min);
	  }
	}
  }

  updateClusterProperties(benlib::range(clusters.size()));

  for (auto& c : clusters) {
	if (c.mass < 1e-5) {
	  std::cout<<"Cluster has mass "<<c.mass<<" and position: "<<c.restCom<<std::endl;
	  return false;
	  exit(0);
	}
  }
  if (!converged) return false;
  std::cout << "numClusters: " << clusters.size() << std::endl;

  return true;
}

