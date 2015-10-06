//#include "accelerationGrid.h"
//template madness

#include <Eigen/Dense>
#include <numeric>
#include <limits>
#include <iostream>

#include "enumerate.hpp"
using benlib::enumerate;
#include "range.hpp"
using benlib::range;


template <typename ParticleType, typename Px>
void AccelerationGrid<ParticleType, Px>::updateGrid(const std::vector<ParticleType>& particles){
  
  assert(numBuckets > 0);
  const auto getter = Px{}
  
  Eigen::Vector3d bbMin = 
	std::accumulate(particles.begin(),
		particles.end(),
		Eigen::Vector3d{std::numeric_limits<double>::max(),
			std::numeric_limits<double>::max(),
			std::numeric_limits<double>::max()},
		[](Eigen::Vector3d best, const ParticleType& p){
		  return Eigen::Vector3d{std::min(best(0),getter(p)(0)),
								 std::min(best(1),getter(p)(1)),
								 std::min(best(2),getter(p)(2))};
		  
		});
  
  
  
  Eigen::Vector3d bbMax = 
	std::accumulate(particles.begin(),
		particles.end(),
		Eigen::Vector3d{std::numeric_limits<double>::min(),
			std::numeric_limits<double>::min(),
			std::numeric_limits<double>::min()},
		[](Eigen::Vector3d best, const ParticleType& p){
		  return Eigen::Vector3d{std::max(best(0),getter(p)(0)),
								 std::max(best(1),getter(p)(1)),
								 std::max(best(2),getter(p)(2))};
		  
		});
  
  delta = (bbMax - bbMin)/numBuckets;
  //give us a little slack space
  bbMin -= 0.05*delta;
  bbMax += 0.05*delta;
  origin = bbMin;
  delta = (bbMax - bbMin)/numBuckets;
  
  grid.assign(numBuckets*numBuckets*numBuckets, {});
  
  for(auto&& e : enumerate(particles)){
	const auto i = e.first;
	const auto& p = e.second;
	const auto bucket = getBucket(getter(p));
	grid[index(bucket)].push_back(i);
  }
}

template <typename ParticleType, typename Px>
std::vector<int> 
AccelerationGrid<ParticleType, Px>::
getNearestNeighbors(const std::vector<ParticleType>& particles, 
					const Eigen::Vector3d& x, double r) const {
  
  //const auto bucket = getBucket(p);
  
  const auto rSquared = r*r;
  //const auto p = Px{}(particles[pIndex]);

  const auto getter = Px{};
  
  std::vector<int> ret;
  
  const auto minBucket = getBucket(x.array() - r);
  const auto maxBucket = getBucket(x.array() + r);
  
  for(auto i : range(minBucket(0), maxBucket(0) +1)){
	for(auto j : range(minBucket(1), maxBucket(1) + 1)){
	  for(auto k : range(minBucket(2), maxBucket(2) +1)){
		//		std::cout << "i " << i << " j " << j << " k " << k << std::endl;
		//		std::cout << "index: " << index(i,j,k) << std::endl;
		//		std::cout << "grid size: " << grid.size() << std::endl;
		for(auto q : grid[index(i,j,k)]){
		  if ((x - getter(particles[q])).squaredNorm() <= rSquared){
			ret.push_back(q);
		  }
		}
	  }
	}
  }
  
  
  return ret;
}
