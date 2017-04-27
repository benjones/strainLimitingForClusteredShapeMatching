#include "world.h"
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

  AccelerationGrid<Particle, EmbeddedPositionGetter> embeddedPositionGrid;
  embeddedPositionGrid.numBuckets = 16; //TODO tune me, 512 buckets now... seems reasonable?
  embeddedPositionGrid.updateGrid(particles);

  
  for(auto& p : particles){p.numClusters = 0; p.totalweight = 0.0;}

  
  while(!lonelyParticles.empty()) {
	std::uniform_int_distribution<> randomDist(0,lonelyParticles.size() -1);
	int r = randomDist(gen);
	//auto currentParticle = lonelyParticles.back();
	//lonelyParticles.pop_back();
	auto currentParticle = lonelyParticles[r];
	lonelyParticles.erase(lonelyParticles.begin() + r);
	Cluster c;
	c.restCom = particles[currentParticle].embeddedPosition;
	std::vector<int> neighbors = embeddedPositionGrid.getNearestNeighbors(particles, c.restCom, neighborRadius);
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
  for (auto& c : clusters) c.cg.init(Eigen::Vector3d::Zero(), neighborRadius);
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


  AccelerationGrid<Particle, EmbeddedPositionGetter> embeddedPositionGrid;
  embeddedPositionGrid.numBuckets = 16; //TODO tune me, 512 buckets now... seems reasonable?
  embeddedPositionGrid.updateGrid(particles);
  
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
	  c.restCom = particles[r].embeddedPosition;
	  clusters.push_back(c);
	}

	converged = false;
	while (!converged) {
	  iters++;
	  converged = true;
	  for (auto& c : clusters) c.members.clear();
	  
	  for (int i=0; i<particles.size(); i++) {
		auto& p = particles[i];
		const auto& pRest = p.embeddedPosition;

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
		embeddedPositionGrid.getNearestNeighbors(particles, c.restCom, params.neighborRadius);

	  c.members.resize(neighbors.size());
	  for(unsigned int i=0; i<neighbors.size(); i++) {
		auto n = neighbors[i];
		Particle &p = particles[n];
		double w = params.kernel(c.restCom - p.embeddedPosition);
		c.members[i] = {n, w};
		p.totalweight += w;
		p.clusters.push_back(j);
		p.numClusters++;
	  }
	}

	for (auto& c : clusters) {
	  c.cg.init(Eigen::Vector3d::Zero(), params.neighborRadius);
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
		  embeddedPositionGrid.getNearestNeighbors(particles, c.restCom, params.neighborRadius);
		c.members.resize(neighbors.size());
		for(unsigned int i=0; i<neighbors.size(); i++) {
		  auto n = neighbors[i];
		  Particle &p = particles[n];
		  double w = params.kernel (c.restCom - p.embeddedPosition);
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
		
		const auto& pRest = p.embeddedPosition;
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
	for (auto& c : clusters) c.cg.init(Eigen::Vector3d::Zero(), params.neighborRadius);
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

  for (auto &c : clusters) {
	for (auto &m : c.members) {
	  m.pos = particles[m.index].embeddedPosition - c.restCom;
	}
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

void World::removeClusters() {
if (clusterFpThreshold < 0) return;
  bool deleteClusters = false;
  // mark cluster for removal
  for (auto &c : clusters) {
	if (c.markedForRemoval || c.newCluster) continue;
	Eigen::JacobiSVD<Eigen::Matrix3d> solver(c.Fp);
	double condFp = solver.singularValues()(0)/solver.singularValues()(2);
	if (condFp > clusterFpThreshold) {
	  c.markedForRemoval = true;
	  c.fadeSteps = clusterFadeOut;
	  c.oweights.reserve(c.members.size());
	  for (auto &m : c.members) c.oweights.push_back(m.weight);
	}
  }

  // fade out dying clusters
  for (auto &c : clusters) {
	if (!c.markedForRemoval || c.fadeSteps <= 0) continue;
	for (auto &&en : benlib::enumerate(c.members)) {
	  auto &m = en.second;
	  double w = c.oweights[en.first]/clusterFadeOut;
	  particles[m.index].totalweight -= w;
	  m.weight -= w;
	}
	c.fadeSteps--;
	if (c.fadeSteps <= 0) deleteClusters = true; 
}  

  if (deleteClusters) {
    std::cout<<"delete Clusters"<<std::endl;
	int sizeBefore = clusters.size();
	std::vector<int> mapping;
	mapping.resize(clusters.size());
	int i1=0, i2=0;
	for (auto &c : clusters) {
	  if (c.markedForRemoval && c.fadeSteps <= 0) {
		mapping[i1++] = -1;
	  } else {
		mapping[i1++] = i2++;
	  }
	}

	for (auto &c : clusters) {
	  std::unordered_set<int> newNeighbors;
	  for (auto &n : c.neighbors) {
		if (mapping[n] != -1)
		  newNeighbors.insert(mapping[n]);
	  }
	  c.neighbors = newNeighbors;
	}
  
	utils::actuallyEraseIf(clusters,
		[](const Cluster& c){
		  return (c.markedForRemoval && c.fadeSteps <= 0);
		}); 
	std::cout << "deleted " << sizeBefore - clusters.size() << " clusters" << std::endl;
  }
  int count = 0;
  for (auto &c : clusters) if (c.markedForRemoval) count++;
  //if (count > 0) std::cout<<count <<" clusters marked for removal"<<std::endl;
}

template <typename Container>
double World::pbdIteration(const Container &particleIndices, const Container &clusterIndices) {
  double omega = 1.0;
  double error = 0.0;
  for(auto& i : particleIndices){
	auto &p = particles[i];
	p.goalPosition.setZero();
	p.tmpd = 0.0;
  }

  for(auto& i : clusterIndices){
	auto &c = clusters[i];
	c.restCom = sumWeightedEmbeddedCOM(c.members);
	
	Eigen::Matrix3d Apq = computeApq(c);
	Eigen::Matrix3d A = Apq * c.aInv * c.Fp.inverse();
	auto pr = utils::polarDecomp(A);
	Eigen::Matrix3d T = pr.first;
	T = T * c.Fp;

	for(const auto& member : c.members){
	  auto &q = particles[member.index];
	  q.goalPosition += member.weight * (T * member.pos + c.restCom);
	  q.tmpd += member.weight;
	}
  }

  for(auto& i : particleIndices){
	auto &p = particles[i];
	p.goalPosition /= p.tmpd;
	error += (p.goalPosition - p.embeddedPosition).squaredNorm();
	p.embeddedPosition = omega*p.goalPosition + (1.0-omega)*p.embeddedPosition;
  }
  return error;
}

void World::addClusters(const ClusteringParams &params) {
  double convergenceThreshold = 1e-6 * params.sqrNeighborRadius; // 0.1% motion allowed

  // loop over particles looking for ones that will not be in any cluster not marked for removal
  // we will seed new clusters at such particles
  for (auto &p : particles) {
	bool unclustered = true;
	for (auto i = p.clusters.begin();  i != p.clusters.end() && unclustered; i++) {
	  if (!clusters[*i].markedForRemoval) unclustered = false;
	}
	if (!unclustered) continue;
	std::cout<<"adding cluster for particle "<<p.id<<std::endl;

	// get a list of particles and clusters to embed in (3D) embedded space
	std::unordered_set<int> particlesToEmbedSet;
	std::unordered_set<int> clustersToEmbedSet;
	for (auto &c : clusters) c.embedId = -1;

	// Try embedding all clusters of p and all neighboring clusters.
	for (auto &i : p.clusters) {
	  auto &c = clusters[i];
	  clustersToEmbedSet.insert(i);
	  for (auto &m : c.members) {
		particlesToEmbedSet.insert(m.index);
		for (auto &j : c.neighbors) {
		  clustersToEmbedSet.insert(j);
		  for (auto &n : clusters[j].members) {
			particlesToEmbedSet.insert(n.index);
		  }
		}
	  }
	}

	std::vector<int> particlesToEmbed(particlesToEmbedSet.begin(), particlesToEmbedSet.end());
	std::vector<int> clustersToEmbed(clustersToEmbedSet.begin(), clustersToEmbedSet.end());
#define LINEARSOLVE 0
#if LINEARSOLVE
	// copy from sets to vectors
	for (unsigned int i = 0; i<clustersToEmbed.size(); i++) clusters[clustersToEmbed[i]].embedId = i;

	// set up the linear system for embedding
	// we will embed particles and clusters
	int psize = particlesToEmbed.size();
	int dim = psize + clustersToEmbed.size();
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(dim,dim);
	Eigen::VectorXd x = Eigen::VectorXd::Zero(dim);
	Eigen::VectorXd y = Eigen::VectorXd::Zero(dim);
	Eigen::VectorXd z = Eigen::VectorXd::Zero(dim);
	Eigen::VectorXd a = Eigen::VectorXd::Zero(dim);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(dim);
	Eigen::VectorXd c = Eigen::VectorXd::Zero(dim);

	// fill in the matrix and rhs
	for (int pindex = 0; pindex<psize; pindex++) {
	  Particle &q = particles[particlesToEmbed[pindex]]; 
	  
	  // loop over all clusters that contain this particle
	  for (auto &j : q.clusters) {
		auto &cluster = clusters[j];
		if (cluster.embedId == -1) continue;
		// look for q in member array
		int k=0;
		while (k < cluster.members.size() && cluster.members[k].index != q.id) {k++;}
		assert (k < cluster.members.size());

		// compute weight
		double w = cluster.members[k].weight/(q.totalweight);
		double w2 = w*w;
		Eigen::Vector3d d = cluster.Fp * cluster.members[k].pos;
		int cindex = psize + cluster.embedId;
		a(pindex) += w*d(0);
		b(pindex) += w*d(1);
		c(pindex) += w*d(2);
		// rhs for cluster is 0.0

		//if (cindex == dim) continue;
		//std::cout<<cindex<<" "<<pindex<<" "<<dim<<" "<<psize<<" "<<cluster.embedId<<std::endl;
		A(pindex, pindex) += w2;
		A(cindex, cindex) += w2;
		A(pindex, cindex) -= w2;
		A(cindex, pindex) -= w2;
	  }
	  if (p.id == q.id) A(pindex,pindex) += 1.0;
	}

	//std::cout<<b<<std::endl;
	  assert(A.allFinite());
	  assert(b.allFinite());

	Eigen::LLT<Eigen::MatrixXd> lltofA(A);
	x = lltofA.solve(a);
	y = lltofA.solve(b);
	z = lltofA.solve(c);

	for (int i = 0; i < psize; i++) {
	  auto &e = particles[particlesToEmbed[i]].embeddedPosition;
	  //std::cout<<"p("<<i<<"): "<<x(i)<<", "<<y(i)<<", "<<z(i)<<std::endl;
	  e(0) = x(i);
	  e(1) = y(i);
	  e(2) = z(i);
	}
	for (int i=0; i<clustersToEmbed.size(); i++) {
	  auto &e = clusters[clustersToEmbed[i]].restCom;
	  //std::cout<<"c("<<i<<"): "<<x(psize+i)<<", "<<y(psize+i)<<", "<<z(psize+i)<<std::endl;
	  e(0) = x(psize+i);
	  e(1) = y(psize+i);
	  e(2) = z(psize+i);
	}
#else
	double initialError = pbdIteration(particlesToEmbedSet, clustersToEmbedSet);
	while (pbdIteration(particlesToEmbedSet, clustersToEmbedSet) > 1e-8*initialError){std::cout<<"."; std::cout.flush();};
	std::cout<<std::endl;
#endif

	//clusters[clustersToEmbed[clustersToEmbed.size()-1]].restCom = Eigen::Vector3d::Zero();
	double embedding_error = 0.0;
	for (auto &i : clustersToEmbedSet) {
	  auto &c = clusters[i];
	  Eigen::Matrix3d Apq = computeApq(c);
	  Eigen::Matrix3d A = Apq * c.aInv * c.Fp.inverse();
	  auto pr = utils::polarDecomp(A);
	  Eigen::Matrix3d T = pr.first;
	  T = T * c.Fp;
	  for (auto &m : c.members) {
		Eigen::Vector3d foo = (T * m.pos) - (particles[m.index].embeddedPosition - c.restCom);
		embedding_error += foo.squaredNorm();
		//std::cout<<foo(0)<<" "<<foo(1)<<" "<<foo(2)<<std::endl;
	  }
	}
	std::cout<<"Embedding Error: "<<embedding_error<<" ("<<clustersToEmbedSet.size()<<" clusters, "<<particlesToEmbedSet.size()<<" particles)"<<std::endl;

	Cluster newCluster;
	newCluster.restCom = particles[p.id].embeddedPosition;
	
	AccelerationGrid<Particle, EmbeddedPositionGetter> embeddedPositionGrid;
	embeddedPositionGrid.numBuckets = 16; //TODO tune me, 512 buckets now... seems reasonable?
	embeddedPositionGrid.updateGrid(particles);
	bool converged = false;

	for (auto &q : particles) {
	  q.futuretotalweight = 0.0;
	}
	for (auto &c : clusters) {
	  if (!c.markedForRemoval) {
		for (auto &m : c.members) {
		  particles[m.index].futuretotalweight += m.weight;
		}
	  }
	}

	while (!converged) {
	  //std::cout<<particles[p.id].embeddedPosition - newCluster.restCom<<std::endl;
	  converged = true;
	  std::vector<int> neighbors =
		embeddedPositionGrid.getNearestNeighbors(particles, newCluster.restCom, params.neighborRadius);
	  newCluster.members.resize(neighbors.size());
	  for(unsigned int i=0; i<neighbors.size(); i++) {
		auto n = neighbors[i];
		Particle &p = particles[n];
		double w = params.kernel (newCluster.restCom - p.embeddedPosition);
		newCluster.members[i] = {n, w};
	  }

	  //compute cluster COMs
	  newCluster.worldCom = newCluster.restCom; //store the last rest COM here for now
	  newCluster.mass = 0.0;
	  newCluster.restCom = Eigen::Vector3d::Zero();
	  for (auto &member : newCluster.members) {
		auto &p = particles[member.index];
		double w = member.weight / (p.futuretotalweight+member.weight);
		newCluster.restCom += w * p.mass * p.embeddedPosition;
		newCluster.mass += w * p.mass;
	  }
	  newCluster.restCom /= newCluster.mass;

	  if ((newCluster.restCom - newCluster.worldCom).squaredNorm() > convergenceThreshold) {
		converged = false;
	  }
	}

	bool found = false;
	newCluster.oweights.resize(newCluster.members.size());
	for (auto &&en : benlib::enumerate(newCluster.members)) {
	  auto &member = en.second;
	  member.pos = particles[member.index].embeddedPosition - newCluster.restCom;
	  newCluster.oweights[en.first] = member.weight;
	  member.weight = 0.0;
	  particles[member.index].clusters.push_back(clusters.size());
	  particles[member.index].numClusters++;
	  if (member.index == p.id) found = true;
	}
	if (!found) std::cout<<"WARNING: Particle not in cluster"<<std::endl;
	//std::cout<<"new cluster: "<<std::endl;
	//for (auto &m : newCluster.members) std::cout<<m.index<<" ";
	//std::cout<<std::endl;
	newCluster.toughness = toughness;
	newCluster.cg.init(Eigen::Vector3d::Zero(), params.neighborRadius);
	newCluster.fadeSteps = clusterFadeIn;
	newCluster.newCluster = true;
	newCluster.Fp = Eigen::Matrix3d::Identity(); 
	newCluster.FpNew = newCluster.Fp;
	brandNewClusters.push_back(clusters.size());
	clusters.push_back(newCluster);
	//std::cout<<"Just set cluster.Fp"<<std::endl<<newCluster.Fp<<std::endl<<newCluster.FpNew<<std::endl;

	for (int i=0; i<clustersToEmbed.size(); i++) {
	  clusters[clustersToEmbed[i]].restCom = Eigen::Vector3d::Zero();
	}
	//std::cout<<A<<std::endl<<computeApq(newCluster)<<std::endl<<computeAqqInv(newCluster)<<std::endl<<std::endl;
  }

  for (auto &c : clusters) {
	if (!c.newCluster || c.fadeSteps <= 0) continue;
	for (auto &&en : benlib::enumerate(c.members)) {
	  auto &m = en.second;
	  double w = c.oweights[en.first]/clusterFadeIn;
	  particles[m.index].totalweight += w;
 	  m.weight += w;
	}
	c.fadeSteps--;
	if (c.fadeSteps <= 0) c.newCluster = false;
	//for (auto &m : c.members) {
	// std::cout<<m.weight<<" "<<m.pos<<std::endl;
	//}
	//exit(-1);
  }  
  int count = 0;
  for (auto &c : clusters) if (c.newCluster) count++;
  //if (count > 0) std::cout<<count <<" new clusters"<<std::endl;
}
