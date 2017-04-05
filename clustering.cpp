#include "clustering.h"
#include "range.hpp"
#include "utils.h"
#include "accelerationGrid.h"

#include "profiler.hpp"
#include <random>


inline double cube(double x) { return x*x*x;}

std::vector<Cluster> makeRandomClusters(std::vector<Particle>& particles, double neighborRadius) {
  // std::random_device rd; unused since we're using the default seed always
  
  std::mt19937 gen(std::mt19937::default_seed);

  
  std::vector<Cluster> clusters;
  
  auto r = benlib::range(particles.size());
  auto lonelyParticles = std::vector<size_t>(r.begin(), r.end());

  AccelerationGrid<Particle, RestPositionGetter> restPositionGrid;
  restPositionGrid.numBuckets = 16; //TODO tune me, 512 buckets now... seems reasonable?
  restPositionGrid.updateGrid(particles);

  
  for(auto& p : particles){p.numClusters = 0; p.totalweight = 0.0;}

  
  while(!lonelyParticles.empty()) {
	std::uniform_int_distribution<> randomDist(0,lonelyParticles.size() -1);
	int r = randomDist(gen);
	//auto currentParticle = lonelyParticles.back();
	//lonelyParticles.pop_back();
	auto currentParticle = lonelyParticles[r];
	lonelyParticles.erase(lonelyParticles.begin() + r);
	Cluster c;
	c.restCom = particles[currentParticle].restPosition;
	std::vector<int> neighbors = restPositionGrid.getNearestNeighbors(particles, c.restCom, neighborRadius);
	c.members.resize(neighbors.size());
	double w = 1.0;
	for(unsigned int i=0; i<neighbors.size(); i++) {
	  auto n = neighbors[i];
	  Particle &p = particles[n];
	  p.numClusters++; 
	  p.totalweight += w; 
	  p.clusters.push_back(clusters.size());
	  c.members[i] = {n, w};
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
	for (auto &member : c.members) {
	  auto &p = particles[member.index];
	  double w = member.weight / p.totalweight;
	  c.restCom += w * p.mass * p.position;
	  c.mass += w * p.mass;
	}
	c.restCom /= c.mass;
  }
  return clusters;
}

std::vector<Cluster> makeClusters(std::vector<Particle>& particles,
	const ClusteringParams& params){
  
  std::cout << "using clusteringAlgorithm "
			<< params.clusteringAlgorithm << std::endl;
  std::vector<Cluster> clusters;

  bool converged;
  int iters = 0;
  double convergenceThreshold = 1e-6 * params.sqrNeighborRadius; // 0.1% motion allowed
  //std::random_device rd; //unused
  std::mt19937 gen(std::mt19937::default_seed);
  std::uniform_int_distribution<> randomDist(0,particles.size() -1);
  std::vector<bool> picked(particles.size(), false);


  if (params.clusteringAlgorithm == 2) {
	return makeRandomClusters(particles, params.neighborRadius);
  }

  AccelerationGrid<Particle, RestPositionGetter> restPositionGrid;
  restPositionGrid.numBuckets = 16; //TODO tune me, 512 buckets now... seems reasonable?
  restPositionGrid.updateGrid(particles);
  
  std::cout << "nClusters: " << params.nClusters << std::endl;
  std::cout << "k means initialization...";
  if (params.clusteringAlgorithm < 2) { //kmeans or fuzzy cmeans

	//seed nClusters clusters
	assert(particles.size() >= params.nClusters);
	while (clusters.size() < params.nClusters) {
	  Cluster c;
	  int r = randomDist(gen);
	  if (picked[r]) continue;
	  picked[r] = true;
	  c.restCom = particles[r].restPosition;
	  clusters.push_back(c);
	}

	converged = false;
	while (!converged) {
	  iters++;
	  converged = true;
	  for (auto& c : clusters) c.members.clear();
	  
	  for (int i=0; i<particles.size(); i++) {
		auto& p = particles[i];
		const auto& pRest = p.restPosition;

		//closest to pRest
		auto bestPair = utils::minProjectedElement(
			clusters,
			[&pRest](const Cluster& c){
			  return (c.restCom - pRest).squaredNorm();
			});
		bestPair.first->members.push_back({i, 1.0});
		int bestCluster = std::distance(clusters.begin(), bestPair.first);
		
		/*
		int bestCluster = 0;
		double bestNorm = (clusters[0].restCom - p.restPosition).squaredNorm();
		for (auto j = 0; j<clusters.size(); j++) {
		  double newNorm = (clusters[j].restCom - p.restPosition).squaredNorm();
		  if (newNorm < bestNorm) {
			bestCluster = j;
			bestNorm = newNorm;
		  }
		  }
		  clusters[bestCluster].members.push_back(std::pair<int,double>(i,1.0));*/
		if (p.numClusters != bestCluster) converged = false;
		p.numClusters = bestCluster;
		p.totalweight = 1.0;
	  }
	  
	  for (auto& c : clusters) {
		c.mass = 0.0;
		c.restCom = Eigen::Vector3d::Zero();
		for (const auto &member : c.members) {
		  auto &p = particles[member.index];
		  double w = member.weight / p.totalweight;
		  c.restCom += w * p.mass * p.position;
		  c.mass += w * p.mass;
		}
		c.restCom /= c.mass;
	  }
	} //end while not converged
  }
  std::cout  << "finished"  << std::endl;
  //fill in weights to clusters and particles
  if (params.clusteringAlgorithm == 1) {
	for(auto& p : particles){p.numClusters = 0; p.totalweight = 0.0;}

	for (auto j=0; j<clusters.size(); j++) {
	  auto &c = clusters[j];
	  std::vector<int> neighbors =
		restPositionGrid.getNearestNeighbors(particles, c.restCom, params.neighborRadius);

	  c.members.resize(neighbors.size());
	  for(unsigned int i=0; i<neighbors.size(); i++) {
		auto n = neighbors[i];
		Particle &p = particles[n];
		double w = params.kernel(c.restCom - p.restPosition);
		c.members[i] = {n, w};
		p.totalweight += w;
		p.clusters.push_back(j);
		p.numClusters++;
	  }
	}

	for (auto& c : clusters) {
	  c.cg.init(c.restCom, params.neighborRadius);
	  c.mass = 0.0;
	  c.restCom = Eigen::Vector3d::Zero();
	  for (auto &member : c.members) {
		auto &p = particles[member.index];
		double w = member.weight / p.totalweight;
		c.restCom += w * p.mass * p.position;
		c.mass += w * p.mass;
	  }
	  c.restCom /= c.mass;
	}
	for (auto& p : particles) {
	  if (p.numClusters == 0) {
		std::cout<<"Particle has no cluster"<<std::endl;
		return {}; //no clusters
	  }
	}
  }

  // fuzzy c-means loop
  if (params.clusteringAlgorithm < 1) {
	std::cout << "starting c means loop..." << std::endl;
	iters = 0;


	while ((!converged || iters < 5) && iters < params.clusterItersMax) {
	  converged = true;
	  iters++;
	  if(iters % 100 == 0){
		std::cout << "starting cmeans iteration " << iters << std::endl;
	  }
	  for (auto& p : particles){
		p.clusters.clear();
		p.numClusters = 0;
		p.totalweight = 0.0;
	  }

	  for (auto j = 0; j<clusters.size(); j++) {
		Cluster &c = clusters[j];
		int oldMembers = c.members.size();
		std::vector<int> neighbors =
		  restPositionGrid.getNearestNeighbors(particles, c.restCom, params.neighborRadius);
		c.members.resize(neighbors.size());
		for(unsigned int i=0; i<neighbors.size(); i++) {
		  auto n = neighbors[i];
		  Particle &p = particles[n];
		  double w = params.kernel (c.restCom - p.restPosition);
		  c.members[i] = {n, w};
		  p.totalweight += w;
		  p.clusters.push_back(j);
		  p.numClusters++;
		}
		if (oldMembers != c.members.size()) {
		  converged = false;
		}
	  }

	  //stick particles not within neighborRadius of a cluster in the nearest cluster
	  for(auto i : range(particles.size())){
		auto& p = particles[i];
		if (p.numClusters > 0) continue;
		
		const auto& pRest = p.restPosition;
		auto bestPair = utils::minProjectedElement(clusters,
			[&pRest](const Cluster& c){ return (c.restCom - pRest).squaredNorm();});
		bestPair.first->members.push_back({static_cast<int>(i), params.blackhole});
		
		/*
		  int bestCluster = 0;
		  double bestNorm = (clusters[0].restCom - p.restPosition).squaredNorm();
		  for (auto j = 0; j<clusters.size(); j++) {
		  double newNorm = (clusters[j].restCom - p.restPosition).squaredNorm();
		  if (newNorm < bestNorm) {
		  bestCluster = j;
		  bestNorm = newNorm;
		  }
		  }
		  clusters[bestCluster].members.push_back(std::pair<int,double>(i,blackhole));
		*/
		p.totalweight = 1.0;
		//std::cout<<"particle "<<i<<" not in cluster"<<std::endl;
		converged = false;
		
	  }

	  //compute the FCM weights x_i_c = 1/(sum_d_inclusters(||x_i - x_c||/||x_i - x_d||)^2/(m-1)
	  if (params.clusterKernel == 4) {

		//almost certain this doesn't need to be nested so deep
		for(const auto i : range(particles.size())){
		  auto& p = particles[i];
		  p.totalweight = 0;
		  //members.second holds distances now
		  std::vector<double> updatedWeights(p.clusters.size());

		  for(auto cIndex : range(p.clusters.size())){
			double sum = 0;
			auto& c = clusters[cIndex];
			auto pInC = utils::findIfInContainer(c.members,
				[i](const Cluster::Member& member){
				  return member.index == i;});
			double xToC = pInC->weight;

			//only the clusters that p is part of
			for(auto& dIndex: p.clusters){
			  auto& d = clusters[dIndex];
			  auto pInD = utils::findIfInContainer(d.members,
				  [i](const Cluster::Member& member){
					return member.index == i;
				  });
			  //it'd better be in there...
			  assert(pInD != d.members.end());
			  sum += pow(xToC/ pInD->weight, 2.0/params.kernelWeight -1);
			}
			updatedWeights[cIndex] = 1.0/sum;
		  }

		  //now overwrite distances in members with weights
		  for(auto&& cIndex : p.clusters){
			auto& c = clusters[cIndex];
			auto pInC = utils::findIfInContainer(c.members,
				[i](const Cluster::Member& member){
				  return member.index == i;
				});
			assert(pInC != c.members.end());
			pInC->weight = updatedWeights[cIndex];
		  }
		  
		}





		for (auto &c : clusters) {
		  for (auto &member : c.members) {
			particles[member.index].totalweight += member.weight;
		  }
		}
	  }


	  //compute cluster COMs
	  for (auto& c : clusters) {
		c.worldCom = c.restCom; //store the last rest COM here for now
		c.mass = 0.0;
		c.restCom = Eigen::Vector3d::Zero();
		for (auto &member : c.members) {
		  auto &p = particles[member.index];
		  double w = member.weight / p.totalweight;
		  c.restCom += w * p.mass * p.position;
		  c.mass += w * p.mass;
		}
		c.restCom /= c.mass;
		if ((c.restCom - c.worldCom).squaredNorm() > convergenceThreshold) {
		  converged = false;
		}
	  }

	}
	for (auto& c : clusters) c.cg.init(c.restCom, params.neighborRadius);
	std::cout<<"Fuzzy c-means clustering ran "<<iters<<" iterations."<<std::endl;
  }

  for (auto& p : particles) {
	if (p.numClusters == 0) {
	  std::cout<<"Particle has no cluster"<<std::endl;
	  return {}; //empty set of clusters means clustering didn't work
	}
  }

  if (!converged) return {}; //no clusters means it didn't work
  std::cout << "numClusters: " << clusters.size() << std::endl;

  for(auto& c : clusters){
    c.initialMembers = c.members.size();
  }

  return clusters;
}


std::vector<Cluster> iterateMakeClusters(
	std::vector<Particle>& particles,
	ClusteringParams& params){

  
  
  
  params.sqrNeighborRadius = params.neighborRadius * params.neighborRadius;
  params.poly6norm = 315.0 / (64.0 * M_PI * cube(cube(params.neighborRadius)));
  std::vector<Cluster> clusters;
  while ((clusters = makeClusters(particles, params)).empty()) {
	params.nClusters = std::ceil(1.25*params.nClusters); 
	params.neighborRadius *= 1.25; 
	if (params.nClusters > params.nClustersMax &&
		params.neighborRadius > params.neighborRadiusMax) {
	  std::cout << "exceeded limits, nclusters: " 
				<< params.nClusters << " " << params.nClustersMax << std::endl
				<< "neighborhood radius: " << params.neighborRadius
				<< " " << params.neighborRadiusMax << std::endl;
	  exit(0);
	}
	params.nClusters = std::min(params.nClusters, params.nClustersMax);
	params.neighborRadius = std::min(params.neighborRadius, params.neighborRadiusMax);
	params.sqrNeighborRadius = params.neighborRadius * params.neighborRadius;
	params.poly6norm = 315.0 / (64.0 * M_PI * cube(cube(params.neighborRadius)));
	std::cout <<"trying clustering again with nClusters = "
			  << params.nClusters <<" neighborRadius = "
			  << params.neighborRadius << std::endl;
  }
  std::cout << "nClusters = " << params.nClusters << std::endl;
  return clusters;

}

  
double ClusteringParams::kernel(const Eigen::Vector3d &x) const{
  switch (clusterKernel) {
  case 0: // 1 / r^2
	return 1.0 / (x.squaredNorm() + 1.0e-4);
  case 1: // constant weight
	return 1.0;
  case 2: // poly 6
	return poly6norm * cube(sqrNeighborRadius-x.squaredNorm());
  case 3: // blend
	return poly6norm * cube(sqrNeighborRadius-x.squaredNorm()) + kernelWeight;
  case 4: // fuzzy c-means
	return x.norm();
  default:
	return 1.0 / (x.squaredNorm() + 1.0e-4);
  }
}
