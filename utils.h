#pragma once

#include <Eigen/Dense>
#include <vector>
#include <utility>

#ifdef __APPLE__
//why, apple?   why????
#include <OpenGL/glu.h>
#else
#include <gl/glu.h>
#endif


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
	
        Eigen::Matrix3d Q = U*V.transpose();
	if(Q.determinant() < 0)
        {
	  sigma(2) *= -1;
	  V.col(2) *= -1;
          Q = U*V.transpose();
	}

        return std::make_pair(Q, V*sigma.asDiagonal()*V.transpose());
  }



  inline void drawSphere(double r, int lats, int longs) {
	int i, j;
	for(i = 1; i <= lats; i++) {
	  double lat0 = M_PI * (-0.5 + (double) (i - 1) / lats);
	  double z0  = sin(lat0);
	  double zr0 =  cos(lat0);
    
	  double lat1 = M_PI * (-0.5 + (double) i / lats);
	  double z1 = sin(lat1);
	  double zr1 = cos(lat1);
    
	  glBegin(GL_QUAD_STRIP);
	  for(j = 0; j <= longs; j++) {
		double lng = 2 * M_PI * (double) j / longs;
		double x = cos(lng);
		double y = sin(lng);
    
		glNormal3f(x * zr1, y * zr1, z1);
		glVertex3f(r*x * zr1, r*y * zr1, r*z1);
		glNormal3f(x * zr0, y * zr0, z0);
		glVertex3f(r*x * zr0, r*y * zr0, r*z0);
	  }
	  glEnd();
	}
  }

  inline void drawClippedSphere(double r, int lats, int longs,
        const Eigen::Vector3d& center,
        const std::vector<std::pair<Eigen::Vector3d, double> >& planes) {
     int i, j;
     for(i = 1; i <= lats; i++) {
        double lat0 = M_PI * (-0.5 + (double) (i - 1) / lats);
        double z0  = sin(lat0);
        double zr0 =  cos(lat0);

        double lat1 = M_PI * (-0.5 + (double) i / lats);
        double z1 = sin(lat1);
        double zr1 = cos(lat1);

        glBegin(GL_QUADS);
        for(j = 0; j < longs; j++) {
           double lng0 = 2 * M_PI * (double) j / longs;
           double x0 = cos(lng0);
           double y0 = sin(lng0);
           double lng1 = 2 * M_PI * ((double) (j+1)) / longs;
           double x1 = cos(lng1);
           double y1 = sin(lng1);

           Eigen::Vector3d n0(x0 * zr1, y0 * zr1, z1);
           Eigen::Vector3d n1(x0 * zr0, y0 * zr0, z0);
           Eigen::Vector3d n2(x1 * zr0, y1 * zr0, z0);
           Eigen::Vector3d n3(x1 * zr1, y1 * zr1, z1);

           Eigen::Vector3d v0 = r*n0+center;
           Eigen::Vector3d v1 = r*n1+center;
           Eigen::Vector3d v2 = r*n2+center;
           Eigen::Vector3d v3 = r*n3+center;

           //always draw if there are no planes
           bool draw = true;

           //if there are planes, check to make sure we're on the right side
           for (int p=0; p<planes.size(); p++) {
              Eigen::Vector3d norm = planes[p].first;
              double offset = planes[p].second;

              if (norm.dot(v0) > -offset &&
                    norm.dot(v1) > -offset &&
                    norm.dot(v2) > -offset &&
                    norm.dot(v3) > -offset) {
                 draw = false;
              }
           }

           if (draw) {
              glNormal3f(n0(0),n0(1),n0(2));
              glVertex3f(v0(0),v0(1),v0(2));
              glNormal3f(n1(0),n1(1),n1(2));
              glVertex3f(v1(0),v1(1),v1(2));
              glNormal3f(n2(0),n2(1),n2(2));
              glVertex3f(v2(0),v2(1),v2(2));
              glNormal3f(n3(0),n3(1),n3(2));
              glVertex3f(v3(0),v3(1),v3(2));
           }
        }
        glEnd();
     }
  }


  inline std::pair<Eigen::Vector3d, Eigen::Vector3d>
	getPlaneTangents(const Eigen::Vector3d& normal){
	
	const Eigen::Vector3d xVector{1,0,0};
	Eigen::Vector3d span1 = (fabs(normal.dot(xVector)) < 0.9) ?
	  normal.cross(xVector) :
	  normal.cross(Eigen::Vector3d{0,1,0});
	Eigen::Vector3d span2 = normal.cross(span1);
	
	span1.normalize();
	span2.normalize();
	
	assert(fabs(normal.squaredNorm() - 1) < 0.01);
	assert(fabs(span1.squaredNorm() - 1) < 0.01);
	assert(fabs(span2.squaredNorm() - 1) < 0.01);
	assert(fabs(normal.dot(span1)) < .01);
	assert(fabs(normal.dot(span2)) < .01);
	assert(fabs(span1.dot(span2)) < .01);
	
	return std::make_pair(span1, span2);
  }



  inline void drawCylinder(const Eigen::Vector3d& supportingPoint,
	  const Eigen::Vector3d& normal,
	  double radius){
	
	const double length = 10;
	const int divisions = 10;
	
	const auto tangents = getPlaneTangents(normal);
	
	glBegin(GL_QUAD_STRIP);
	for(int i = 0; i <= divisions; ++i){
	  double theta = i*2.0*M_PI/divisions;
	  Eigen::Vector3d offset = radius*(std::cos(theta)*tangents.first + std::sin(theta)*tangents.second);
	  Eigen::Vector3d v1 = supportingPoint - length*normal + offset;
	  Eigen::Vector3d v2 = supportingPoint + length*normal + offset;
	  glVertex3dv(v1.data());
	  glVertex3dv(v2.data());
	}
	glEnd();

  }

  template <typename Cont, typename Value>
  void actuallyErase(Cont& container, const Value& value){
	container.erase(std::remove(container.begin(), container.end(), value),
		container.end());
  }

  template <typename Cont, typename Pred>
  void actuallyEraseIf(Cont& container, Pred&& predicate){
	container.erase(
		std::remove_if(container.begin(), container.end(), std::forward<Pred>(predicate)),
		container.end());
  }

  inline void drawPlane(const Eigen::Vector3d& normal, double offset, double size) {
     Eigen::Vector3d tangent1, tangent2;

     tangent1 = normal.cross(Eigen::Vector3d{1,0,0});
     if(tangent1.isZero(1e-3)){
        tangent1 = normal.cross(Eigen::Vector3d{0,0,1});
        if(tangent1.isZero(1e-3)){
           tangent1 = normal.cross(Eigen::Vector3d{0,1,0});
        }
     }
     tangent1.normalize();

     tangent2 = normal.cross(tangent1);
     tangent2.normalize(); //probably not necessary

     const double sos = normal.dot(normal);
     const Eigen::Vector3d supportPoint{normal.x()*offset/sos,
        normal.y()*offset/sos,
        normal.z()*offset/sos};

     glBegin(GL_QUADS);
     glNormal3dv(normal.data());
     glVertex3dv((supportPoint + size*(tangent1 + tangent2)).eval().data());
     glVertex3dv((supportPoint + size*(-tangent1 + tangent2)).eval().data());
     glVertex3dv((supportPoint + size*(-tangent1 - tangent2)).eval().data());
     glVertex3dv((supportPoint + size*(tangent1  - tangent2)).eval().data());
     glEnd();
  }

  template <typename Cont, typename Value>
  bool containsValue(const Cont& container, const Value& value){
	return std::find(container.begin(), container.end(), value) != container.end();
  }
  
  inline void drawPlane(const Eigen::Vector3d& normal,
	  double offset, double size, const Eigen::Vector3d &center) {
	Eigen::Vector3d tangent1, tangent2;
	
	tangent1 = normal.cross(Eigen::Vector3d{1,0,0});
     if(tangent1.isZero(1e-3)){
	   tangent1 = normal.cross(Eigen::Vector3d{0,0,1});
        if(tangent1.isZero(1e-3)){
           tangent1 = normal.cross(Eigen::Vector3d{0,1,0});
        }
     }
     tangent1.normalize();

     tangent2 = normal.cross(tangent1);
     tangent2.normalize(); //probably not necessary

     //const double sos = normal.dot(normal);
	 const double d = normal.dot(center) + offset;
	 const Eigen::Vector3d supportPoint = center - d * normal;
     //const Eigen::Vector3d supportPoint{normal.x()*offset/sos,
	 //  normal.y()*offset/sos,
	 //  normal.z()*offset/sos};

     glBegin(GL_QUADS);
     glNormal3dv(normal.data());
     glVertex3dv((supportPoint + size*(tangent1 + tangent2)).eval().data());
     glVertex3dv((supportPoint + size*(-tangent1 + tangent2)).eval().data());
     glVertex3dv((supportPoint + size*(-tangent1 - tangent2)).eval().data());
     glVertex3dv((supportPoint + size*(tangent1  - tangent2)).eval().data());
     glEnd();
  }


  //pick the element of the container that minimizes Functor(elem), for example, the distance to a point
  template <typename Cont, typename Functor,
			typename It = typename Cont::iterator,
			typename Val = decltype(std::declval<Functor>()(*std::declval<It>()))>
  std::pair<It, Val>
  minProjectedElement(Cont& cont, Functor&& functor){
	using std::begin;
	using std::end;
	if(begin(cont) == end(cont)){
	  return std::make_pair(end(cont), std::numeric_limits<Val>::quiet_NaN());
	}

	It bestIt = begin(cont);
	Val bestVal = functor(*bestIt);

	for(auto it = begin(cont); it != end(cont); ++it){
	  Val thisVal = functor(*it);
	  if(thisVal < bestVal){
		bestVal = thisVal;
		bestIt = it;
	  }
	}
	return std::make_pair(bestIt, bestVal);
  }

  template<typename Cont, typename Val,
		   typename It = typename Cont::iterator>
  It findInContainer(Cont& cont, const Val& val){
	using std::begin;
	using std::end;
	return std::find(begin(cont), end(cont), val);
  }
  template<typename Cont, typename Pred,
		   typename It = typename Cont::iterator>
  It findIfInContainer(Cont& cont, Pred&& pred){
	using std::begin;
	using std::end;
	return std::find_if(begin(cont), end(cont), std::forward<Pred>(pred));
  }

  template<typename T>
  struct incomplete;
}
