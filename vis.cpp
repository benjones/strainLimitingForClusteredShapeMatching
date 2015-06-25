#include "world.h"
#include <fstream>

#ifdef __APPLE__
//why, apple?   why????
#include <OpenGL/glu.h>
#else
#include <gl/glu.h>
#endif

#include "utils.h"
#include "color_spaces.h"

#include "enumerate.hpp"
using benlib::enumerate;
#include "range.hpp"
using benlib::range;


void World::draw(SDL_Window* window) const {
  throw("don't call me");
  glEnable(GL_DEPTH_TEST);
  glFrontFace(GL_CCW);
  glEnable(GL_CULL_FACE);
  glClearColor(0.2, 0.2, 0.2, 1);
  glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  int windowWidth, windowHeight;
  SDL_GetWindowSize(window, &windowWidth, &windowHeight);
  gluPerspective(45,static_cast<double>(windowWidth)/windowHeight,
				 .5, 100);

  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(cameraPosition.x(), cameraPosition.y(), cameraPosition.z(),
			cameraLookAt.x(), cameraLookAt.y(), cameraLookAt.z(),
			cameraUp.x(), cameraUp.y(), cameraUp.z());



  drawPlanes();
  
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glDisable(GL_DEPTH_TEST);
  //draw clusters
  glMatrixMode(GL_MODELVIEW);
  if(drawClusters){
	for(auto&& pr : benlib::enumerate(clusters)){
	  auto& c = pr.second;
	  auto i = pr.first;
	  glPushMatrix();

	  auto com = sumWeightedWorldCOM(c.neighbors);
	  glTranslated(com.x(), com.y(), com.z());
	  glColor4d(i/(2.0*clusters.size()), 1.0, 1.0, 0.3);
	  utils::drawSphere(c.renderWidth, 10, 10);
	  glPopMatrix();
	}
  }
  glEnable(GL_DEPTH_TEST);


  
  if(!particles.empty()){						
	//	glDisable(GL_DEPTH_TEST);
	glColor3f(1,1,1);
	glPointSize(5);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE,
					sizeof(Particle),
					&(particles[0].position));
	glDrawArrays(GL_POINTS, 0, particles.size());
	/*	glPointSize(3);
	glColor3f(0,0,1);
	glVertexPointer(3, GL_DOUBLE, sizeof(Particle),
					&(particles[0].goalPosition));
	glDrawArrays(GL_POINTS, 0, particles.size());
	*/

  }


  glFlush();
  SDL_GL_SwapWindow(window);
}



void World::drawPretty(SDL_Window* window) const {

  glEnable(GL_DEPTH_TEST);
  glFrontFace(GL_CCW);
  glEnable(GL_CULL_FACE);
  glClearColor(0.2, 0.2, 0.2, 1);
  glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  int windowWidth, windowHeight;
  SDL_GetWindowSize(window, &windowWidth, &windowHeight);
  gluPerspective(45,static_cast<double>(windowWidth)/windowHeight,
				 .5, 100);

  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(cameraPosition.x(), cameraPosition.y(), cameraPosition.z(),
			cameraLookAt.x(), cameraLookAt.y(), cameraLookAt.z(),
			cameraUp.x(), cameraUp.y(), cameraUp.z());



  drawPlanesPretty();
  
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


  auto max_t = 0.0;
  if (colorByToughness) {
     for(auto&& pr : benlib::enumerate(clusters)){
        auto& c = pr.second;
        if (c.toughness < std::numeric_limits<double>::infinity()) {
           if (c.toughness > max_t) {
              max_t = c.toughness;
           }
        }
     }
  }


  glDisable(GL_DEPTH_TEST);
  //draw clusters
  glMatrixMode(GL_MODELVIEW);
  if(drawClusters){
     for(auto&& pr : benlib::enumerate(clusters)){
        auto& c = pr.second;
        const auto i = pr.first;
        if (which_cluster == -1 || i == which_cluster) {
           glPushMatrix();

           auto com = sumWeightedWorldCOM(c.neighbors);
           glTranslated(com.x(), com.y(), com.z());
           RGBColor rgb = HSLColor(2.0*acos(-1)*(i%12)/12.0, 0.7, 0.7).to_rgb();
           //RGBColor rgb = HSLColor(2.0*acos(-1)*i/clusters.size(), 0.7, 0.7).to_rgb();
           if (colorByToughness) {
              if (c.toughness == std::numeric_limits<double>::infinity()) {
                 rgb = RGBColor(0.0, 0.0, 0.0);
              } else {
                 auto factor = c.toughness/max_t;
                 rgb = RGBColor(1.0-factor, factor, factor);
              }
           }
           glColor4d(rgb.r, rgb.g, rgb.b, 0.3);
           utils::drawSphere(c.renderWidth, 10, 10);
           glPopMatrix();
        }
     }
  }
  glEnable(GL_DEPTH_TEST);


  glDisable(GL_DEPTH_TEST);
  //draw fracture planes 
  glMatrixMode(GL_MODELVIEW);
  if(drawFracturePlanes){
     for(auto&& pr : benlib::enumerate(clusters)){
        auto& c = pr.second;
        const auto i = pr.first;
        if (which_cluster == -1 || i == which_cluster) {
           glPushMatrix();

           auto& cg = c.cg;

           Eigen::Matrix4d gl_trans = Eigen::Matrix4d::Identity();
           gl_trans.block<3,3>(0,0) << c.restToWorldTransform;
           if (cg.planes.size() > 0) {
              std::cout << gl_trans << std::endl;
           }

           RGBColor rgb = HSLColor(2.0*acos(-1)*(i%12)/12.0, 0.7, 0.7).to_rgb();
           glColor4d(rgb.r, rgb.g, rgb.b, 0.3);
           glPopMatrix();
        }
     }
  }
  glEnable(GL_DEPTH_TEST);

  if (which_cluster != -1) {
     auto& c = clusters[which_cluster];

	 Eigen::Matrix3d Apq = computeApq(c);
	 Eigen::Matrix3d A = Apq*c.aInv;
	 if (nu > 0.0) A = A*c.Fp.inverse(); // plasticity
	 
	 //do the SVD here so we can handle fracture stuff
	 Eigen::JacobiSVD<Eigen::Matrix3d> solver(A, 
		 Eigen::ComputeFullU | Eigen::ComputeFullV);
	 
	 Eigen::Matrix3d U = solver.matrixU(), V = solver.matrixV();
	 Eigen::Vector3d sigma = solver.singularValues();

	 std::cout << "sigma: " << sigma << std::endl;

     glColor4d(0,0,0, 0.9);

     glPointSize(5);

     if (!drawColoredParticles) {
        glBegin(GL_POINTS);
        for(auto &i : c.neighbors){
           glVertex3dv(particles[i.first].position.data());
        }
        glEnd();
     }
  }

  
  if(!particles.empty()){						
     //	glDisable(GL_DEPTH_TEST);
     glPointSize(10);

     if (!drawColoredParticles) {
        glColor4d(1,1,1,0.8);
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_DOUBLE,
              sizeof(Particle),
              &(particles[0].position));
        glDrawArrays(GL_POINTS, 0, particles.size());
     } else {
        glBegin(GL_POINTS);
        for(auto& p : particles) {
           //find nearest cluster

           int min_cluster = p.clusters[0];
           auto com = clusters[min_cluster].worldCom;//computeNeighborhoodCOM(clusters[min_cluster]);
           Eigen::Vector3d dir = p.position - com;
           double min_sqdist = dir.squaredNorm();
           for (auto& cInd : p.clusters) {
			 com = clusters[cInd].worldCom;//computeNeighborhoodCOM(clusters[cInd]);
              dir = p.position - com;
              double dist = dir.squaredNorm();
              if (dist < min_sqdist) {
                 min_sqdist = dist;
                 min_cluster = cInd;
              }
           }

           RGBColor rgb = HSLColor(2.0*acos(-1)*(min_cluster%12)/12.0, 0.7, 0.7).to_rgb();
           //RGBColor rgb = HSLColor(2.0*acos(-1)*min_cluster/clusters.size(), 0.7, 0.7).to_rgb();
           if ((which_cluster == -1 || min_cluster == which_cluster) &&
                 clusters[min_cluster].neighbors.size() > 1) {
           //      sqrt(min_sqdist) < 0.55*clusters[min_cluster].renderWidth) {
              glColor4d(rgb.r, rgb.g, rgb.b, 0.8);
           } else {
              glColor4d(1,1,1,0.8);
           }
           //glPushMatrix();
           //glTranslated(p.position[0], p.position[1], p.position[2]);
           //utils::drawSphere(0.01, 4, 4);
           glVertex3dv(p.position.data());
           //glPopMatrix();
        }
        glEnd();
     }


     /*	glPointSize(3);
         glColor3f(0,0,1);
         glVertexPointer(3, GL_DOUBLE, sizeof(Particle),
         &(particles[0].goalPosition));
         glDrawArrays(GL_POINTS, 0, particles.size());
      */

  }


  glColor3d(1, 0, 1);
  for(auto& projectile : projectiles){
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	auto currentPosition = projectile.start + elapsedTime*projectile.velocity;
	glTranslated(currentPosition.x(), currentPosition.y(), currentPosition.z());
	utils::drawSphere(projectile.radius, 10, 10);
	glPopMatrix();
  }

  glColor3d(0,1,0);
  for(auto& cylinder : cylinders){
	utils::drawCylinder(cylinder.supportPoint, cylinder.normal, cylinder.radius);
  }

  glFlush();
  SDL_GL_SwapWindow(window);
}


void World::drawSingleCluster(SDL_Window* window, int frame) const {

  glEnable(GL_DEPTH_TEST);
  glFrontFace(GL_CCW);
  glEnable(GL_CULL_FACE);
  glClearColor(0.2, 0.2, 0.2, 1);
  glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  int windowWidth, windowHeight;
  SDL_GetWindowSize(window, &windowWidth, &windowHeight);
  gluPerspective(45,static_cast<double>(windowWidth)/windowHeight,
				 .5, 100);

  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(cameraPosition.x(), cameraPosition.y(), cameraPosition.z(),
			cameraLookAt.x(), cameraLookAt.y(), cameraLookAt.z(),
			cameraUp.x(), cameraUp.y(), cameraUp.z());



  drawPlanes();

  glDisable(GL_DEPTH_TEST);
  //draw clusters
  glMatrixMode(GL_MODELVIEW);
  auto i = frame % clusters.size();
  auto& c = clusters[i];

  glPushMatrix();

  auto com = sumWeightedWorldCOM(c.neighbors);
  glTranslated(com.x(), com.y(), com.z());
  glColor4d(static_cast<double>(i)/clusters.size(), 0, 1.0, 0.3);
  utils::drawSphere(c.width, 10, 10);
  glPopMatrix();
  glEnable(GL_DEPTH_TEST);


  
  //	glDisable(GL_DEPTH_TEST);
  glColor3f(1,1,1);
  glPointSize(3);
  
  glBegin(GL_POINTS);
  for(auto &i : c.neighbors){
	glVertex3dv(particles[i.first].position.data());
  }
  glEnd();

  glFlush();
  SDL_GL_SwapWindow(window);
}


void World::drawPlanes() const{
  //draw planes
  glDisable(GL_CULL_FACE);
  double totalCount = planes.size() + movingPlanes.size() + twistingPlanes.size() + tiltingPlanes.size();
  for(auto&& pr : enumerate(planes)){
	const auto i = pr.first;
	const auto& plane = pr.second;
	glColor4d(0.5, static_cast<double>(i)/totalCount,
			  0.5, 1);

	drawPlane(plane.head(3), plane.w());
  }
  glDepthMask(false);
  for(auto&& pr : enumerate(movingPlanes)){
	const auto i = pr.first + planes.size();
	const auto& plane = pr.second;
	glColor4d(0.5, i/totalCount, 0.5, 1);
	drawPlane(plane.normal, plane.offset + elapsedTime*plane.velocity);
  }
  for(auto&& pr : enumerate(twistingPlanes)){
	const auto i = pr.first + planes.size() + movingPlanes.size();
	const auto& plane = pr.second;
	glColor4d(0.5, i/totalCount, 0.5, 1);
	if (elapsedTime <= plane.lifetime)
	  drawTPlane(plane.normal, plane.offset, elapsedTime*plane.angularVelocity, plane.width);
  }
  for(auto&& pr : enumerate(tiltingPlanes)){
	const auto i = pr.first + planes.size() + movingPlanes.size() + twistingPlanes.size();
	const auto& plane = pr.second;
	glColor4d(0.5, i/totalCount, 0.5, 1);
	if (elapsedTime <= plane.lifetime)
	  drawTiltPlane(plane.normal, plane.tilt, plane.offset, elapsedTime*plane.angularVelocity, plane.width);
  }

  glDepthMask(true);
}


void World::drawPlanesPretty() const{
  //draw planes
  glDisable(GL_CULL_FACE);
  for(auto&& pr : enumerate(planes)){
	const auto i = pr.first;
	const auto& plane = pr.second;
   RGBColor rgb = HSLColor(0.25*acos(-1)*i/planes.size()+0.0*acos(-1), 0.3, 0.7).to_rgb();
	glColor4d(rgb.r, rgb.g, rgb.b, 1.0);

	drawPlane(plane.head(3), plane.w());
  }
  for(auto&& pr : enumerate(movingPlanes)){
	const auto i = pr.first; 
	const auto& plane = pr.second;
   RGBColor rgb = HSLColor(0.25*acos(-1)*i/movingPlanes.size()+1.0*acos(-1), 0.3, 0.7).to_rgb();
	glColor4d(rgb.r, rgb.g, rgb.b, 1.0);
	drawPlane(plane.normal, plane.offset + elapsedTime*plane.velocity);
  }
  for(auto&& pr : enumerate(twistingPlanes)){
	const auto i = pr.first; 
	const auto& plane = pr.second;
   RGBColor rgb = HSLColor(0.25*acos(-1)*i/twistingPlanes.size()+0.5*acos(-1), 0.3, 0.7).to_rgb();
	glColor4d(rgb.r, rgb.g, rgb.b, 1.0);
	if (elapsedTime <= plane.lifetime)
	  drawTPlane(plane.normal, plane.offset, elapsedTime*plane.angularVelocity, plane.width);
  }
  for(auto&& pr : enumerate(tiltingPlanes)){
	const auto i = pr.first; 
	const auto& plane = pr.second;
   RGBColor rgb = HSLColor(0.25*acos(-1)*i/tiltingPlanes.size()+1.5*acos(-1), 0.3, 0.7).to_rgb();
	glColor4d(rgb.r, rgb.g, rgb.b, 1.0);
	if (elapsedTime <= plane.lifetime)
	  drawTiltPlane(plane.normal, plane.tilt, plane.offset, elapsedTime*plane.angularVelocity, plane.width);
  }



}


void World::drawPlane(const Eigen::Vector3d& normal, double offset) const{
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


	
	const double size = 100;
	glBegin(GL_QUADS);
	glNormal3dv(normal.data());
	glVertex3dv((supportPoint + size*(tangent1 + tangent2)).eval().data());
	glVertex3dv((supportPoint + size*(-tangent1 + tangent2)).eval().data());
	glVertex3dv((supportPoint + size*(-tangent1 - tangent2)).eval().data());
	glVertex3dv((supportPoint + size*(tangent1  - tangent2)).eval().data());
	glEnd();

}

void World::drawTPlane(const Eigen::Vector3d& normal, double offset, double roffset, double width) const{
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
	
   Eigen::AngleAxisd t(roffset,normal);
   tangent1 = t * tangent1; 
   tangent2 = t * tangent2; 

	const double sos = normal.dot(normal);
	const Eigen::Vector3d supportPoint{normal.x()*offset/sos,
		normal.y()*offset/sos,
		normal.z()*offset/sos};


	
	const double size = width;
	glBegin(GL_QUADS);
	glNormal3dv(normal.data());
	glVertex3dv((supportPoint + size*(tangent1 + tangent2)).eval().data());
	glVertex3dv((supportPoint + size*(-tangent1 + tangent2)).eval().data());
	glVertex3dv((supportPoint + size*(-tangent1 - tangent2)).eval().data());
	glVertex3dv((supportPoint + size*(tangent1  - tangent2)).eval().data());
	glEnd();

}


void World::drawTiltPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& tilt, double offset, double roffset, double width) const{
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
	
   Eigen::AngleAxisd t(roffset,tilt);
   tangent1 = t * tangent1; 
   tangent2 = t * tangent2; 

	const double sos = normal.dot(normal);
	const Eigen::Vector3d supportPoint{normal.x()*offset/sos,
		normal.y()*offset/sos,
		normal.z()*offset/sos};


	
	const double size = width;
	glBegin(GL_QUADS);
	glNormal3dv(normal.data());
	glVertex3dv((supportPoint + size*(tangent1 + tangent2)).eval().data());
	glVertex3dv((supportPoint + size*(-tangent1 + tangent2)).eval().data());
	glVertex3dv((supportPoint + size*(-tangent1 - tangent2)).eval().data());
	glVertex3dv((supportPoint + size*(tangent1  - tangent2)).eval().data());
	glEnd();

}

void World::zoom(int amount){
  if(amount < -1){
	cameraPosition -= 0.1*(cameraLookAt - cameraPosition);
  } else if(amount > 1){
	cameraPosition += 0.1*(cameraLookAt - cameraPosition);
  }
}

void World::pan(Eigen::Vector2i oldPosition, Eigen::Vector2i newPosition){

  const Eigen::Vector2d delta = (newPosition - oldPosition).eval().template cast<double>();
  
  const Eigen::Vector3d forwardVector = cameraLookAt - cameraPosition;
  const Eigen::Vector3d rightVector = forwardVector.cross(cameraUp);
  const Eigen::Vector3d upVector = rightVector.cross(forwardVector);

  const double scale = 0.0005;
  
  cameraPosition += scale*(-delta.x()*rightVector +
						   delta.y()*upVector);


}

void World::move(bool forward){
  const double scale = 0.01* (forward ? 1 : -1);
  
  const Eigen::Vector3d delta = scale*(cameraLookAt - cameraPosition);
  
  cameraLookAt += delta;
  cameraPosition += delta;
  
}



