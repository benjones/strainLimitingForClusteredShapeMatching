#include "world.h"
#include <fstream>

#include "utils.h"

#include <random>
#include <iostream>

#include "accelerationGrid.h"

#include "enumerate.hpp"
using benlib::enumerate;
#include "range.hpp"
using benlib::range;

inline double sqr (const double &x) {return x*x;}

void World::timestep(){

  for(auto& twistingPlane : twistingPlanes){
	for(auto& p : particles){
      if (twistingPlane.outside(p) && twistingPlane.lifetime < elapsedTime) 
		p.outsideSomeMovingPlane = false;
   }
  }

  for(auto& tiltingPlane : tiltingPlanes){
	for(auto& p : particles){
      if (tiltingPlane.outside(p) && tiltingPlane.lifetime < elapsedTime) 
		p.outsideSomeMovingPlane = false;
   }
  }


  //scope block for profiler
  std::vector<FractureInfo> potentialSplits;
  {
	auto timer = prof.timeName("dynamics");
	for(auto& p : particles){
	  p.oldPosition = p.position;
	  p.goalPosition.setZero();
	  p.goalVelocity.setZero();
	}
	
	
	
	for(auto&& en : benlib::enumerate(clusters)){
	  auto& cluster = en.second;
	  cluster.worldCom = sumWeightedWorldCOM(cluster.members);
	  Eigen::Vector3d clusterVelocity = computeClusterVelocity(cluster);
	  
	  Eigen::Matrix3d Apq = computeApq(cluster);
	  Eigen::Matrix3d A = Apq*cluster.aInv;
	  if (nu > 0.0) A = A*cluster.Fp.inverse(); // plasticity
	  
	  //do the SVD here so we can handle fracture stuff
	  Eigen::JacobiSVD<Eigen::Matrix3d> solver(A, 
		  Eigen::ComputeFullU | Eigen::ComputeFullV);
	  
	  Eigen::Matrix3d U = solver.matrixU(), V = solver.matrixV();
	  Eigen::Vector3d sigma = solver.singularValues();
	  
	  //std::cout << "sigma " << sigma << std::endl;

	  if(fabs(sigma(0) - 1.0) > cluster.toughness){
		potentialSplits.emplace_back(en.first, 
			sigma(0) - cluster.toughness, 
			V.col(0));
	  }

	  {
		auto timer = prof.timeName("plasticity");
		cluster.updatePlasticity(sigma, U, V, yield, nu, hardening);
	  }

	  Eigen::Matrix3d T = U*V.transpose();
	  if (nu > 0.0) T = T*cluster.Fp; // plasticity
	  
	  for(auto &member : cluster.members){
		auto &p = particles[member.index];
		double w = member.weight / p.totalweight;
		p.goalPosition += w*(T*(p.restPosition - cluster.restCom) + cluster.worldCom);
		p.goalVelocity += w*clusterVelocity;
	  }	  
	}
	
	
	assertFinite();
	
	for(auto& p : particles){
	  if(p.numClusters > 0){
		//p.goalPosition /= p.numClusters;
		//p.goalVelocity /= p.numClusters;
	  } else {
		p.goalPosition = p.position;
		p.goalVelocity = p.velocity;
	  }
	  if (!p.outsideSomeMovingPlane) {
		p.velocity += dt * gravity + (alpha/dt)*(p.goalPosition- p.position) + 
		  springDamping*(p.goalVelocity - p.velocity); 
	  }
	  p.position += dt * p.velocity;
	}
	
	assertFinite();

	for(auto iter : range(numConstraintIters)){
	  (void)iter; //unused
	  strainLimitingIteration();
	}

	for(auto& p : particles){
	  p.velocity = (1.0/dt)*(p.position - p.oldPosition);
	}
	for (auto &c : clusters) c.worldCom = sumWeightedWorldCOM(c.members);	
	assertFinite();
  }

  //std::cout<<"doFracture"<<std::endl;
  if (fractureOn) {
	doFracture(std::move(potentialSplits));
	//std::cout<<"splitoutliers"<<std::endl;
	splitOutliers();
	//std::cout<<"cullsmallclusters"<<std::endl;
	cullSmallClusters();
	//std::cout<<"removelonelyparticles"<<std::endl;
	removeLonelyParticles();
	//std::cout<<"updateClusterProperties"<<std::endl;
  }

  updateClusterProperties(range(clusters.size()));

  //std::cout<<"updateTransforms"<<std::endl;
  for (auto &c : clusters) updateTransforms(c);
  //std::cout<<"selfCollisions"<<std::endl;
  if (selfCollisionsOn) selfCollisions();

  bounceOutOfPlanes();

  for(auto& c : clusters){
	c.Fp = c.FpNew;
	c.worldCom = sumWeightedWorldCOM(c.members);	
	c.renderWidth = 0;
	for(auto& member : c.members){
	  c.renderWidth = std::max(c.renderWidth, (c.worldCom - particles[member.index].position).norm());
	}
  }

  elapsedTime += dt;
  //std::cout << "elapsed time: " << elapsedTime << std::endl;
}

void World::strainLimitingIteration(){
  const double gammaSquared = gamma * gamma; 
  for(auto& p : particles){
		p.goalPosition.setZero();
  }
  for(auto& c : clusters){
	c.worldCom = sumWeightedWorldCOM(c.members);
	
	Eigen::Matrix3d Apq = computeApq(c);
	Eigen::Matrix3d A = Apq*c.aInv;
	if (nu > 0.0) A = A*c.Fp.inverse(); // plasticity
	auto pr = utils::polarDecomp(A);

	Eigen::Matrix3d T = pr.first;
	if (nu > 0.0) T = T * c.Fp;

	for(const auto& member : c.members){
	  auto &q = particles[member.index];
	  Eigen::Vector3d rest = (q.restPosition - c.restCom);
	  Eigen::Vector3d goal = T*(rest) + c.worldCom;
	  double ratio = (goal-q.position).squaredNorm() / (c.width*c.width);
	  
	  if (ratio > gammaSquared) {
		q.goalPosition += 
		  (member.weight/q.totalweight) * (goal + 
			  sqrt(gammaSquared/ratio) * 
			  (q.position - goal));
	  } else {
		q.goalPosition += (member.weight/q.totalweight) * q.position;
	  }
	}
	
  }
  
  for(auto& p : particles){
	if(p.numClusters == 0) {
	  p.goalPosition = p.position;
	}
	p.position = omega*p.goalPosition + (1.0-omega)*p.position;
  }
}

///////////////////////////////
// FRACTURE
///////////////////////////////

void World::doFracture(std::vector<World::FractureInfo> potentialSplits){
  auto timer = prof.timeName("fracture");
  //do fracture
  std::sort(potentialSplits.begin(), potentialSplits.end(),
	  [](const FractureInfo& a, const FractureInfo& b){
		return std::get<1>(a) < std::get<1>(b);
	  });

  if(!potentialSplits.empty()){
	std::cout << "potential splits: " << potentialSplits.size() << std::endl;
  }

  bool aCancel = false;
  for(auto &ps : potentialSplits){
	size_t cIndex = std::get<0>(ps);

	auto worldCOM = clusters[cIndex].worldCom;

	//recompute A matrix:
	Eigen::Matrix3d Apq = computeApq(clusters[cIndex]);
	Eigen::Matrix3d A = Apq*clusters[cIndex].aInv;
	if (nu > 0.0) A = A*clusters[cIndex].Fp.inverse(); // plasticity
	
	//do the SVD here so we can handle fracture stuff
	Eigen::JacobiSVD<Eigen::Matrix3d> solver(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
	
	Eigen::Vector3d sigma = solver.singularValues();
	if(fabs(sigma(0) - 1) < clusters[cIndex].toughness){
	  if(!aCancel){
		aCancel = true;
		std::cout << "cancelled fracture with updated stuff" << std::endl;
	  }
	  continue;
	}
	
	Eigen::Vector3d splitDirection = solver.matrixV().col(0);

	auto& c = clusters[cIndex];	  

	auto it = std::partition(c.members.begin(), c.members.end(),
		[&worldCOM, &splitDirection, this](const Cluster::Member& a){
		  return (worldCOM - particles[a.index].position).dot(splitDirection) > 0;
		});


	auto oldSize = std::distance(c.members.begin(), it);
	auto newSize = std::distance(it, c.members.end());
	if(newSize == 0 || oldSize == 0){ continue;}

	//expected to be mostly in teh x direction for the cube example, and it was
	//std::cout << "split direction: " << splitDirection << std::endl;
	
	//make a new cluster
	Cluster newCluster;
	newCluster.members.resize(newSize);
	for(auto i : range(newSize)){
	  newCluster.members[i] = *(it + i);
	}

	// copy relevant variables
	newCluster.Fp = c.Fp; // plasticity
	newCluster.FpNew = c.FpNew;
	newCluster.cstrain = c.cstrain; // plasticity
	newCluster.toughness = c.toughness;

	//delete the particles from the old one
	c.members.erase(c.members.begin() + oldSize, c.members.end());

	// Update Collision Geometry
	newCluster.cg = c.cg;
	//Eigen::Matrix3d T = solver.matrixU()*solver.matrixV().transpose();
	//if (nu > 0.0) T = T*c.Fp; // plasticity
	//T = T.inverse().eval();
	//Eigen::Vector3d n = T*splitDirection;
	//we should be transforming from world to rest, so we shouldn't include plasticity
	Eigen::Vector3d n = ((Apq*clusters[cIndex].aInv).inverse()*splitDirection).normalized();

	//we need to worry about signs at some point
	c.cg.addPlane(n, -(n.dot(c.restCom)));
	newCluster.cg.addPlane(-n, (n.dot(c.restCom)));
	
	clusters.push_back(newCluster);	  

	updateClusterProperties(std::initializer_list<size_t>{cIndex, clusters.size()-1});
	
	{
	  auto timerTwo = prof.timeName("propagate");
	  //split from other clusters
	  std::vector<int> allParticles;
	  for (auto &member : clusters[cIndex].members) allParticles.push_back(member.index);
	  for (auto &member : newCluster.members) allParticles.push_back(member.index);

	  std::vector<size_t> affectedClusters; //keep sorted

	  // Adam found a segfault from an invalidated iterator here on 6/29 and applied a "quick fix" since Ben was going to redo this logic anyway.
	  std::vector<Particle> newParticles;
	  for(auto& member : allParticles){
		auto& particle = particles[member];
		for(auto thisIndex : particle.clusters){
		  //insert into sorted array
		  auto it = std::lower_bound(affectedClusters.begin(), affectedClusters.end(), thisIndex);
		  if(it == affectedClusters.end() || *it != thisIndex){
			affectedClusters.insert(it, thisIndex);
		  }
		  
		  auto& thisCluster = clusters[thisIndex];

		  if (thisIndex >= clusters.size()) {
			std::cout<<particle.clusters.size()<<std::endl;
			std::cout<<thisIndex<<">="<<clusters.size()<<std::endl<<member<<" "<<particles.size()<<std::endl;
			for (auto foo : particle.clusters) {
			  std::cout<<foo<<" ";
			}
			std::cout<<std::endl;
		  }
		  assert(thisIndex < clusters.size());

		  if(((particle.position - worldCOM).dot(splitDirection) >= 0) !=
			  ((thisCluster.worldCom - worldCOM).dot(splitDirection) >= 0 )){
			unsigned int i = 0;
			while (thisCluster.members[i].index != member && i < thisCluster.members.size()) i++;
			if (i == thisCluster.members.size()) break;

			double w = thisCluster.members[i].weight;

			Particle q(particle);
			q.clusters.clear();
			q.clusters.push_back(thisIndex);
			q.numClusters = 1;
			q.mass = (w / particle.totalweight) * particle.mass;
			q.totalweight = w;
			double newMass = ((particle.totalweight - w) / particle.totalweight) * particle.mass;
			//if (newMass < 0.1*particle.mass || q.mass < 0.1*particle.mass) {
			  // in this case we should just delete the particle from the cluster and let the mass be lost...
			  //std::cout<<"mass low "<<n<<" "<<newMass<<" "<<q.mass<<" "<<p.numClusters<<std::endl;
			  //continue;
			//}

			// This seems like a good idea, but it invalidates the iterator to particle.clusters
			//particle.clusters.erase(
			//	std::remove(particle.clusters.begin(), particle.clusters.end(),
			//		thisIndex), particle.clusters.end());
			particle.numClusters--;
			particle.mass = newMass;
			particle.totalweight -= w;
			particle.flags |= Particle::SPLIT;
			q.flags |= Particle::SPLIT;

			newParticles.push_back(q);
			thisCluster.members[i].index = particles.size()+newParticles.size()-1;
		  }
		}
	  }
	  particles.insert(particles.end(),newParticles.begin(), newParticles.end());
	  updateClusterProperties(affectedClusters);
	
	  //	for(auto& c : affectedClusters){
	  //	  clusters[c].toughness *= 0.995;
	//	}
	
	//break;
	}
	
  }
}	

void World::splitOutliers() {
  for(auto&& en : benlib::enumerate(clusters)){
	auto& c = en.second;
	bool updateCluster = false;
	for(auto &member : c.members) {
	  auto &p = particles[member.index];
	  if ((p.numClusters > 1) && ((p.position - c.worldCom).norm() > 
			  (1.0 + gamma) * outlierThreshold * (p.restPosition - c.restCom).norm())) {
		Particle q(p);
		q.clusters.clear();
		q.clusters.push_back(en.first);
		q.numClusters = 1;
		q.mass = (member.weight/p.totalweight) * p.mass;
		q.totalweight = member.weight;
		double newMass = ((p.totalweight-member.weight)/p.totalweight) * p.mass;
		//if (newMass < 0.1*p.mass || q.mass < 0.1*p.mass) {
		//continue;
		//}

		utils::actuallyErase(p.clusters, en.first);

		p.numClusters--;
		p.mass = newMass;
		p.totalweight -= member.weight;

		p.flags |= Particle::SPLIT;
		q.flags |= Particle::SPLIT;

		particles.push_back(q);
		member.index = particles.size()-1;
		updateCluster = true;
	  }
	} 
  }
}

void World::cullSmallClusters() {
  double threshold = 1e-4;
  auto sizeBefore = clusters.size();
  for (auto &c : clusters) {
	if (c.members.size() < 4 || c.mass < threshold) {
	  for (auto &member : c.members) {
		particles[member.index].totalweight -= member.weight;
	  }
	}
  }


  utils::actuallyEraseIf(clusters,
	  [](const Cluster& c){
		return (c.members.size() < 4 || c.mass < 1e-4);
	  }); 

  if(clusters.size() != sizeBefore){
	std::cout << "deleted " << sizeBefore - clusters.size() << " clusters" << std::endl;
  }
  countClusters();
}

void World::removeLonelyParticles() {
  std::vector<int> mapping;
  mapping.resize(particles.size());
  int i1=0, i2=0;
  for (auto &p : particles) {
	if (p.numClusters > 0) {
	  mapping[i1++] = i2++;
	} else {
	  mapping[i1++] = -1;
	}
  }

  for (auto &c : clusters) {
	for (auto &member : c.members) {
	  assert (mapping[member.index] != -1);
	  member.index = mapping[member.index];
	}
  }
  
  utils::actuallyEraseIf(particles,
	  [](const Particle& p){
		return (p.numClusters == 0);
	  });
}

/////////////////////////////////////
// COLLISIONS
/////////////////////////////////////

void World::bounceOutOfPlanes(){
  prof.timeName("ground collisions");
  bool bounced = true;
  int iters = 0;
  
  const double epsilon = 1e-5;
  
  while(bounced){
	bounced = false;
	
	for(auto & plane : planes){
	  Eigen::Vector3d norm = plane.head(3);
	  for(auto & p : particles){
		if(p.position.dot(norm) < plane.w()){
		  bounced = true;
		  p.position += (epsilon + plane.w() - p.position.dot(norm))*norm;
		  
		  //zero velocity in the normal direction
		  p.velocity -= p.velocity.dot(norm)*norm;
		  p.velocity *= 0.4; //friction
		  
		}
		
	  }
	}
	
	++iters;
  }

   //handle moving planes
  for(auto& movingPlane : movingPlanes){
	for(auto& p : particles){
      if (dragWithPlanes) {
   	  movingPlane.dragParticle(p, elapsedTime);
      } else {
   	  movingPlane.bounceParticle(p, elapsedTime);
      }
   }
  }

   //handle twisting planes
  for(auto& twistingPlane : twistingPlanes){
	for(auto& p : particles){
   	  twistingPlane.twistParticle(p, elapsedTime);
   }
  }

   //handle tilting planes
 for(auto& tiltingPlane : tiltingPlanes){
	for(auto& p : particles){
   	  tiltingPlane.tiltParticle(p, elapsedTime);
   }
  }


  for(auto& projectile : projectiles){
	for(auto& particle : particles){
	  projectile.bounceParticle(particle, elapsedTime);
	}
  }

  for(auto& cylinder : cylinders){
	for(auto& particle : particles){
	  cylinder.bounceParticle(particle);
	}
  }
}

void World::buildClusterMaps() {
  countClusters();
  std::vector<std::vector<int > > idToClusters;
  idToClusters.resize(particles.size());
  for (auto &p : particles) {
	idToClusters[p.id].insert(idToClusters[p.id].end(), p.clusters.begin(), p.clusters.end());
  }

  clusterCollisionMap.clear();
  clusterCollisionMap.resize(clusters.size());

  for (auto &i : idToClusters) {
	for (auto &c : i) {
	  for (auto &d : i) {
		if (c == d) continue;
		if (!utils::containsValue(clusterCollisionMap[c], d)){
		  clusterCollisionMap[c].push_back(d);
		}
	  }
	}
  }
	//for (auto &p : particles) {
	//for (auto &c : p.clusters) {
	//	  for (auto &d : p.clusters) {
	//		if (c == d) continue;
	//		if (!utils::containsValue(clusterCollisionMap[c], d)){
	//		  clusterCollisionMap[c].push_back(d);
	//		}
	//	  }
	//	}
	// }
}

bool CollisionGeometry::project(Eigen::Vector3d &x) {
  //return false;
  Eigen::Vector3d y;
  double n;
  Eigen::Vector3d d = x - c;
  double m = d.norm();
  if (m >= r) return false;

  y = c + (r/m)*d;
  n = r-m;

  for (auto &p : planes) {
	m = p.first.dot(x) + p.second;
  	if (m >= 0.0) return false;
  	if (-m < n) {
  	  n = -m;
  	  y = x - m*p.first;
  	}
  }

  assert((x-y).norm() <= n+1e-4);
  x = y;
  return true;
}

inline Eigen::Vector3d restToWorld(const Cluster &c, const Eigen::Vector3d &x) { 
  assert(x.allFinite());
  Eigen::Vector3d y = c.worldCom + c.restToWorldTransform * (x-c.restCom);
  assert(y.allFinite());
  return y;
}


inline Eigen::Vector3d worldToRest(const Cluster &c, const Eigen::Vector3d &x) {
  assert(x.allFinite());
  Eigen::Vector3d y = c.restCom + c.worldToRestTransform * (x-c.worldCom);
  if (!y.allFinite()) {
	std::cout<<x<<std::endl;
	std::cout<<y<<std::endl;
	std::cout<<c.worldToRestTransform<<std::endl;
	std::cout<<c.members.size()<<std::endl;
  }
  assert(y.allFinite());
  return y;
}

struct ClusterComGetter{
  Eigen::Vector3d operator()(const Cluster& c) const{
	return c.worldCom;
  }
};

struct ClusterRadiusGetter{
  double operator()(const Cluster& c) const{
	return c.width;
  }
};


void World::selfCollisions() {
  const double alpha = collisionRestitution;
  buildClusterMaps();

  AccelerationGrid<Cluster, ClusterComGetter> accelerationGrid;
  accelerationGrid.numBuckets = 16; //todo, tune me
  accelerationGrid.updateGridWithRadii(clusters, ClusterRadiusGetter{});
  auto potentialClusterPairs = accelerationGrid.getPotentialPairs();
  
  for(auto& clusterPair : potentialClusterPairs){
	auto i = clusterPair.first;
	auto j = clusterPair.second;

	if( utils::containsValue(clusterCollisionMap[i], j)){ continue; }

	auto& c = clusters[i];
	auto& d = clusters[j];
	
	if (d.members.size() < 10) continue;
	if ((c.worldCom - d.worldCom).squaredNorm() < sqr(c.width + d.width)) {

	  for (auto& member : c.members){
		auto &p = particles[member.index];
		//if (p.flags & Particle::SPLIT) continue;
		if ((p.position - d.worldCom).squaredNorm() < sqr(d.width)) {
		  // these next two lines look unnecessary with clustermaps
		  // it = find(p.clusters.begin(), p.clusters.end(), j);
		  // (it == p.clusters.end()) {
		  if(!utils::containsValue(p.clusters, j)){ //easier to read... --Ben
			Eigen::Vector3d x = worldToRest(d, p.position);
			if (d.cg.project(x)) {
			  x = restToWorld(d, x);

			  Eigen::Vector3d y = d.worldCom + d.width * (p.position - d.worldCom).normalized();

			  if ((x-p.position).squaredNorm() > (y-p.position).squaredNorm()) {x=y;}
			  
			  p.position = alpha*x + (1.0-alpha)*p.position;
			}
		  }
		}
	  }

	  //swap c and d
	  for (auto& member : d.members){
		auto &p = particles[member.index];
		//if (p.flags & Particle::SPLIT) continue;
		if ((p.position - c.worldCom).squaredNorm() < sqr(c.width)) {
		  // these next two lines look unnecessary with clustermaps
		  if(!utils::containsValue(p.clusters, i)){
			Eigen::Vector3d x = worldToRest(c, p.position);
			if (c.cg.project(x)) {
			  x = restToWorld(c, x);
			  Eigen::Vector3d y = c.worldCom + c.width * (p.position - c.worldCom).normalized();
			  if ((x-p.position).squaredNorm() > (y-p.position).squaredNorm()) {x=y;}
			  
			  p.position = alpha*x + (1.0-alpha)*p.position;
			  
			}
		  }
		}
	  }
	}
  }
	/*
  for (auto && en1 : benlib::enumerate(clusters)) {
	auto &c = en1.second;
	if (c.members.size() < 10) continue;
	for (auto && en2 : benlib::enumerate(clusters)) {
	  if (en1.first == en2.first) continue;
	  // 6/24: I think the old == was a bug here
	  // 6/29, this helper function should be easier to use than raw std::find --Ben
	  if( utils::containsValue(clusterCollisionMap[en1.first], en2.first)){ continue; }

	  auto &d = en2.second;
	  if (d.members.size() < 10) continue;
	  if ((c.worldCom - d.worldCom).squaredNorm() < sqr(c.width + d.width)) {
		for (auto& member : c.members){
		  auto &p = particles[member.index];
		  //if (p.flags & Particle::SPLIT) continue;
		  if ((p.position - d.worldCom).squaredNorm() < sqr(d.width)) {
			// these next two lines look unnecessary with clustermaps
			auto it = find(p.clusters.begin(), p.clusters.end(), en2.first);
			if (it == p.clusters.end()) {
			  Eigen::Vector3d x = worldToRest(d, p.position);
			  if (d.cg.project(x)) {
				//std::cout<<en1.first<<" "<<en2.first<<" "<<n.first<<std::endl;
				//for (auto foo : clusterCollisionMap[en1.first]) std::cout<<foo<<" ";
				//std::cout<<std::endl;
				//for (auto foo : clusterCollisionMap[en2.first]) std::cout<<foo<<" ";
				//std::cout<<std::endl;
				x = restToWorld(d, x);
				Eigen::Vector3d y = d.worldCom + d.width * (p.position - d.worldCom).normalized();
				if ((x-p.position).squaredNorm() > (y-p.position).squaredNorm()) {x=y;}
				//if ((x-p.position).norm() > (p.position-d.worldCom).norm()) {
				if ((x-p.position).norm() > d.width) {
				  std::cout<<"--------------------------"
						   <<(x-p.position).norm()<<" "<<d.width<<" "<<d.cg.r<<std::endl;
				  std::cout<<d.worldCom(0)<<" "<<d.worldCom(1)<<" "<<d.worldCom(2)<<std::endl;
				  std::cout<<p.position(0)<<" "<<p.position(1)<<" "<<p.position(2)<<std::endl;
				  std::cout<<x(0)<<" "<<x(1)<<" "<<x(2)<<std::endl;
				}
				//std::cout<<"before: "<<(p.position - d.worldCom).norm()<<" "<<d.width<<std::endl;
				p.position = alpha*x + (1.0-alpha)*p.position;
				//std::cout<<"after: "<<(p.position - d.worldCom).norm()<<" "<<d.width<<std::endl;
			  }
			}
		  }
		}
	  }
	}
	}*/
}






///////////////////////////////////////////
// UTILS
///////////////////////////////////////////

Eigen::Matrix3d World::computeApq(const Cluster& c) const{
  Eigen::Matrix3d Apq;
  Apq.setZero();
  for (auto &member : c.members) {
	auto &p = particles[member.index];
	Eigen::Vector3d pj = p.position - c.worldCom;
	Eigen::Vector3d qj = p.restPosition - c.restCom;
	Apq += (member.weight/p.totalweight)*p.mass * pj * qj.transpose();
  }
  return Apq;
} 	

Eigen::Vector3d World::computeClusterVelocity(const Cluster &c) const {
  double mass = 0.0;
  Eigen::Vector3d vel = Eigen::Vector3d::Zero();
  for (auto &member : c.members) {
	auto &p = particles[member.index];
	double w = (member.weight / p.totalweight) * p.mass;
	mass += w;
	vel += w * p.velocity;
  }
  return (vel / mass);
}

void World::updateTransforms(Cluster& c) const{
  Eigen::Matrix3d A = computeApq(c)*c.aInv;
  //if (nu > 0.0) A = A*c.Fp.inverse(); 

  //auto pr = utils::polarDecomp(A);
  //Eigen::Matrix3d T = pr.first;
  //if (nu > 0.0) T = T*c.Fp; 
  c.restToWorldTransform = A;
  c.worldToRestTransform = A.inverse();
} 	

void World::countClusters(){
  for(auto& p : particles){
	p.numClusters = 0;
	p.clusters.clear();
  }
  for(auto cInd : benlib::range(clusters.size())){
	for(auto& member : clusters[cInd].members){
	  auto &p = particles[member.index];
	  ++(p.numClusters);
	  p.clusters.push_back(cInd);
	}
  }
  for(auto& p : particles){assert(p.numClusters == p.clusters.size());}
}

// updates coms, width, mass, aInv
template <typename Container>
void World::updateClusterProperties(const Container& clusterIndices){
  countClusters(); //TODO, just touch a subset...

  for(auto& p : particles){
	if (p.mass <= 0.0) std::cout<<p.mass<<" "<<std::endl;
	assert(p.mass > 0.0);
  }

  // compute cluster mass, com, width, and aInv
  for(auto cIndex : clusterIndices){
	auto& c = clusters[cIndex];
	//c.Fp.setIdentity(); // plasticity
	c.mass = sumWeightedMass(c.members);
	if (!(c.mass >= 0)) {
	  std::cout<<c.mass<<" "<<c.members.size()<<" "<<c.members[0].weight<<std::endl;
	}
	assert(c.mass >= 0);
	c.restCom = sumWeightedRestCOM(c.members);
	c.worldCom = sumWeightedWorldCOM(c.members);

	if(!c.restCom.allFinite()){std::cout << c.restCom << std::endl;}
	if(!c.worldCom.allFinite()){std::cout << c.worldCom << std::endl;}

	assert(c.restCom.allFinite());
	assert(c.worldCom.allFinite());
	c.width = 0.0;
	for(const auto& member : c.members){
	  c.width = std::max(c.width, (c.restCom - particles[member.index].restPosition).norm());
	} 
	c.renderWidth = c.width; //it'll get updated soon enough
	//assert(c.width >= 0);
	
	c.aInv.setZero();  
	// adam says: this should take weights into account
	for (const auto& member : c.members) {
	  auto &p = particles[member.index];
	  Eigen::Vector3d qj = p.restPosition - c.restCom;
	  c.aInv += (member.weight/p.totalweight)*p.mass * qj * qj.transpose();
	}
	  
	
	//do pseudoinverse
	Eigen::JacobiSVD<Eigen::Matrix3d> solver(c.aInv, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Vector3d sigInv;
	for(auto i : range(3)){
	  if (solver.singularValues()(i) < 1e-6) c.Fp.setIdentity();
	  //if (solver.singularValues()(i) < 1e-6) std::cout<<"yikes "<<solver.singularValues()(i)<<std::endl;
	  sigInv(i) = fabs(solver.singularValues()(i)) > 1e-12 ? 1.0/solver.singularValues()(i) : 0;
	}
	c.aInv = solver.matrixV()*sigInv.asDiagonal()*solver.matrixU().transpose();//c.aInv.inverse().eval();
	if(!c.aInv.allFinite()){
	  std::cout << c.aInv << std::endl;
	  std::cout << solver.singularValues() << std::endl;
	}
	assert(c.aInv.allFinite());
  }

}


