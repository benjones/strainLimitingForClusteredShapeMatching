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



  inline void drawSphere(double r, int lats, int longs) {
	int i, j;
	for(i = 0; i <= lats; i++) {
	  double lat0 = M_PI * (-0.5 + (double) (i - 1) / lats);
	  double z0  = sin(lat0);
	  double zr0 =  cos(lat0);
    
	  double lat1 = M_PI * (-0.5 + (double) i / lats);
	  double z1 = sin(lat1);
	  double zr1 = cos(lat1);
    
	  glBegin(GL_QUAD_STRIP);
	  for(j = 0; j <= longs; j++) {
		double lng = 2 * M_PI * (double) (j - 1) / longs;
		double x = cos(lng);
		double y = sin(lng);
    
		glNormal3f(x * zr0, y * zr0, z0);
		glVertex3f(r*x * zr0, r*y * zr0, r*z0);
		glNormal3f(x * zr1, y * zr1, z1);
		glVertex3f(r*x * zr1, r*y * zr1, r*z1);
	  }
	  glEnd();
	}
  }

  
}
