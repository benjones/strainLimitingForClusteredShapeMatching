#pragma once

#include <Eigen/Dense>
#include <utility>

namespace utils{

  //R, S
  /*
	Replacing this with something faster is very likely to give you a big speed boost
   */
  inline std::pair<Eigen::Matrix3d, Eigen::Matrix3d> polarDecomp(const Eigen::Matrix3d& A){
	
	Eigen::JacobiSVD<Eigen::Matrix3d> solver(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
	
	Eigen::Matrix3d U = solver.matrixU();
	Eigen::Matrix3d V = solver.matrixV();
	Eigen::Vector3d sigma = solver.singularValues();
	
	if(sigma.prod() < 0){
	  sigma(2) *= -1;
	  V.col(2) *= -1;
	}

	return std::make_pair(U*V.transpose(), V*sigma.asDiagonal()*V.transpose());
  }
  
}
