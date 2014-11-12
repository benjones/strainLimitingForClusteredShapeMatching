#pragma once

#include <vector>
#include <Eigen/Dense>


//Px(P) returns an Eigen::Vector3d with the position of P
template <typename ParticleType, typename Px>
class AccelerationGrid{

public:
  std::vector<std::vector<int>> grid;
  
  std::vector<int> getNearestNeighbors(const std::vector<ParticleType>& particles, 
									   const Eigen::Vector3d& x, double radius) const;
  
  void updateGrid(const std::vector<ParticleType>& particles);
  
  Eigen::Vector3d origin, delta;
  int numBuckets;
  
  Eigen::Vector3i inline getBucket(Eigen::Vector3d position) const{
	return Eigen::Vector3i {((position - 
							  origin).array()/
							 delta.array()).template cast<int>()
		}.unaryExpr([this](int a){
			return std::max(0,std::min(a, numBuckets -1));
		  });
  }
  
  
  inline int index(int i, int j, int k) const {
	return i*numBuckets*numBuckets +
	  j*numBuckets + k;
  }
  
  inline int index(Eigen::Vector3i v) const {
	return index(v(0), v(1), v(2));
  }
  
};


//template madness
#include "accelerationGrid.cpp"
