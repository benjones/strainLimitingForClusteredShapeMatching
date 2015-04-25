#include "world.h"
#include <fstream>

#ifdef __APPLE__
//why, apple?   why????
#include <OpenGL/glu.h>
#else
#include <gl/glu.h>
#endif

#include "json/json.h"
#include "utils.h"
#include "color_spaces.h"

#include <random>

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

	  auto com = computeNeighborhoodCOM(c);
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

           auto com = computeNeighborhoodCOM(c);
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



  if (which_cluster != -1) {
     auto& c = clusters[which_cluster];

	 Eigen::Matrix3d init;
	 init.setZero();
	 
	 
	 Eigen::Matrix3d Apq = computeApq(c, init, c.worldCom);
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
        for(auto i : c.neighbors){
           glVertex3dv(particles[i].position.data());
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

  auto com = computeNeighborhoodCOM(c);
  glTranslated(com.x(), com.y(), com.z());
  glColor4d(static_cast<double>(i)/clusters.size(), 0, 1.0, 0.3);
  utils::drawSphere(c.width, 10, 10);
  glPopMatrix();
  glEnable(GL_DEPTH_TEST);


  
  //	glDisable(GL_DEPTH_TEST);
  glColor3f(1,1,1);
  glPointSize(3);
  
  glBegin(GL_POINTS);
  for(auto i : c.neighbors){
	glVertex3dv(particles[i].position.data());
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



void World::loadFromJson(const std::string& _filename){
  elapsedTime = 0;
  filename = _filename;

  std::ifstream ins(filename);
  
  Json::Value root;
  Json::Reader jReader;

  if(!jReader.parse(ins, root)){
	std::cout << "couldn't read input file: " << filename << '\n'
			  << jReader.getFormattedErrorMessages() << std::endl;
	exit(1);
  }


  auto particleFilesIn = root["particleFiles"];
  for(auto i : range(particleFilesIn.size())){
	readParticleFile(particleFilesIn[i].asString());
  }
  
  auto cameraPositionIn = root["cameraPosition"];
  if(!cameraPositionIn.isNull() && cameraPositionIn.isArray() && cameraPositionIn.size() == 3){
	cameraPosition.x() = cameraPositionIn[0].asDouble();
	cameraPosition.y() = cameraPositionIn[1].asDouble();
	cameraPosition.z() = cameraPositionIn[2].asDouble();
  } else {
	cameraPosition = Eigen::Vector3d{0,7, -5};
  }

  auto cameraLookAtIn = root["cameraLookAt"];
  if(!cameraLookAtIn.isNull() && cameraLookAtIn.isArray() && cameraLookAtIn.size() == 3){
	cameraLookAt.x() = cameraLookAtIn[0].asDouble();
	cameraLookAt.y() = cameraLookAtIn[1].asDouble();
	cameraLookAt.z() = cameraLookAtIn[2].asDouble();
  } else {
	cameraLookAt = Eigen::Vector3d::Zero();
  }
  
  auto cameraUpIn = root["cameraUp"];
  if(!cameraUpIn.isNull() && cameraUpIn.isArray() && cameraUpIn.size() == 3){
	cameraUp.x() = cameraUpIn[0].asDouble();
	cameraUp.y() = cameraUpIn[1].asDouble();
	cameraUp.z() = cameraUpIn[2].asDouble();
  } else {
	cameraUp = Eigen::Vector3d{0,1,0};
  }

  dt = root.get("dt",1/60.0).asDouble();
  neighborRadius = root.get("neighborRadius", 0.1).asDouble();
  nClusters = root.get("nClusters", -1).asInt();
  numConstraintIters = root.get("numConstraintIters", 5).asInt();
  alpha = root.get("alpha", 1.0).asDouble();
  omega = root.get("omega", 1.0).asDouble();
  gamma = root.get("maxStretch", 0.1).asDouble();
  gamma = root.get("gamma", gamma).asDouble();
  springDamping = root.get("springDamping", 0.0).asDouble();
  toughness = root.get("toughness", std::numeric_limits<double>::infinity()).asDouble();

  // plasticity parameters
  yield = root.get("yield", 0.0).asDouble();
  nu = root.get("nu", 0.0).asDouble();
  hardening = root.get("hardening", 0.0).asDouble();
  


  auto& gravityIn = root["gravity"];
  if(!gravityIn.isNull() && gravityIn.isArray() && gravityIn.size() == 3){
	gravity.x() = gravityIn[0].asDouble();
	gravity.y() = gravityIn[1].asDouble();
	gravity.z() = gravityIn[2].asDouble();
  } else {
	std::cout << "default gravity" << std::endl;
	gravity = Eigen::Vector3d{0, -9.81, 0};
  }

  auto& planesIn =  root["planes"];
  for(auto i : range(planesIn.size())){
	if(planesIn[i].size() != 4){
	  std::cout << "not a good plane... skipping" << std::endl;
	  continue;
	}

	//x, y, z, (of normal), then offset value
	planes.push_back(Eigen::Vector4d{planesIn[i][0].asDouble(),
		  planesIn[i][1].asDouble(),
		  planesIn[i][2].asDouble(),
		  planesIn[i][3].asDouble()});
	planes.back().head(3).normalize();
  }


  auto& movingPlanesIn = root["movingPlanes"];
  for(auto i : range(movingPlanesIn.size())){
	auto normalIn = movingPlanesIn[i]["normal"];
	if(normalIn.size() != 3){
	  std::cout << "bad moving plane, skipping" << std::endl;
	  continue;
	}
	Eigen::Vector3d normal(normalIn[0].asDouble(),
		normalIn[1].asDouble(),
		normalIn[2].asDouble());
   normal.normalize();
	movingPlanes.emplace_back(normal, 
		movingPlanesIn[i]["offset"].asDouble(),
		movingPlanesIn[i]["velocity"].asDouble());

  }
  std::cout << movingPlanes.size() << " moving planes" << std::endl;

  auto& twistingPlanesIn = root["twistingPlanes"];
  for(auto i : range(twistingPlanesIn.size())){
	auto normalIn = twistingPlanesIn[i]["normal"];
	if(normalIn.size() != 3){
	  std::cout << "bad twisting plane, skipping" << std::endl;
	  continue;
	}
	Eigen::Vector3d normal(normalIn[0].asDouble(),
		normalIn[1].asDouble(),
		normalIn[2].asDouble());
   normal.normalize();
	twistingPlanes.emplace_back(normal, 
		twistingPlanesIn[i]["offset"].asDouble(),
		twistingPlanesIn[i]["angularVelocity"].asDouble(),
		twistingPlanesIn[i]["width"].asDouble(),
		twistingPlanesIn[i].get("lifetime", std::numeric_limits<double>::max()).asDouble());

  }
  std::cout << twistingPlanes.size() << " twisting planes" << std::endl;

  auto& tiltingPlanesIn = root["tiltingPlanes"];
  for(auto i : range(tiltingPlanesIn.size())){
	auto normalIn = tiltingPlanesIn[i]["normal"];
	if(normalIn.size() != 3){
	  std::cout << "bad tilting plane, skipping" << std::endl;
	  continue;
	}
	Eigen::Vector3d normal(normalIn[0].asDouble(),
		normalIn[1].asDouble(),
		normalIn[2].asDouble());
   normal.normalize();
   auto tiltIn = tiltingPlanesIn[i]["tilt"];
	if(tiltIn.size() != 3){
	  std::cout << "bad tilting plane, skipping" << std::endl;
	  continue;
	}
	Eigen::Vector3d tilt(tiltIn[0].asDouble(),
		tiltIn[1].asDouble(),
		tiltIn[2].asDouble());
   tilt.normalize();

	tiltingPlanes.emplace_back(normal, tilt,
		tiltingPlanesIn[i]["offset"].asDouble(),
		tiltingPlanesIn[i]["angularVelocity"].asDouble(),
		tiltingPlanesIn[i]["width"].asDouble(),
		tiltingPlanesIn[i].get("lifetime", std::numeric_limits<double>::max()).asDouble());

  }
  std::cout << tiltingPlanes.size() << " tilting planes" << std::endl;

  auto& projectilesIn = root["projectiles"];
  for(auto i : range(projectilesIn.size())){
	auto& projectile = projectilesIn[i];
	Eigen::Vector3d start(projectile["start"][0].asDouble(),
		projectile["start"][1].asDouble(),
		projectile["start"][2].asDouble());
	Eigen::Vector3d vel(projectile["velocity"][0].asDouble(),
		projectile["velocity"][1].asDouble(),
		projectile["velocity"][2].asDouble());

	projectiles.emplace_back(start, vel, projectile["radius"].asDouble(),
		projectile.get("momentumScale", 0).asDouble());

  }
  
  auto& cylindersIn = root["cylinders"];
  for(auto i : range(cylindersIn.size())){
	auto& cylinder = cylindersIn[i];
	Eigen::Vector3d normal(cylinder["normal"][0].asDouble(),
		cylinder["normal"][1].asDouble(),
		cylinder["normal"][2].asDouble());

	Eigen::Vector3d supportPoint(cylinder["supportPoint"][0].asDouble(),
		cylinder["supportPoint"][1].asDouble(),
		cylinder["supportPoint"][2].asDouble());


	cylinders.emplace_back(normal, supportPoint, cylinder["radius"].asDouble());
  }


  double mass = root.get("mass", 0.1).asDouble();
  for(auto& p : particles){ p.mass = mass;}

  for(auto& p : particles){ p.outsideSomeMovingPlane = false;}

  for(auto& movingPlane : movingPlanes){
	for(auto& p : particles){
      p.outsideSomeMovingPlane |= movingPlane.outside(p);
   }
  }

  for(auto& twistingPlane : twistingPlanes){
	for(auto& p : particles){
      p.outsideSomeMovingPlane |= twistingPlane.outside(p);
   }
  }

  for(auto& tiltingPlane : tiltingPlanes){
	for(auto& p : particles){
      p.outsideSomeMovingPlane |= tiltingPlane.outside(p);
   }
  }



  
  restPositionGrid.numBuckets = 6;
  restPositionGrid.updateGrid(particles);
  
  makeClusters();

  //apply initial rotation/scaling, etc
  //todo, read this from json
  //Eigen::Vector3d axis{0,0,1};
  //for (auto &p : particles) {
	//p.position.x() *= 2.5;
	//p.velocity = 0.5*p.position.cross(axis);
	//auto pos = p.position;
	//.position[0] = 0.36*pos[0] + 0.48*pos[1] - 0.8*pos[2];
	//p.position[1] = -0.8*pos[0] + 0.6*pos[1];// - 0.8*pos[2];
	//p.position[2] = 0.48*pos[0] + 0.64*pos[1] + 0.6*pos[2];
  //}
}




void World::readParticleFile(const std::string& _filename){


  std::ifstream ins(_filename);
  
  Eigen::Vector3d pos;
  ins >> pos.x() >> pos.y() >> pos.z();
  //auto xpos = pos;
  //pos[0] = pos[0]/2.0;
  //pos[1] = pos[1]/2.0 + 2;
  //pos[2] = pos[2]/2.0;
  double bbMin[3] = {pos.x(), pos.y(), pos.z()}, 
	bbMax[3] = {pos.x(), pos.y(), pos.z()};
  while(ins){
	particles.emplace_back();
	particles.back().position = pos;
	particles.back().restPosition = pos;
	particles.back().velocity = Eigen::Vector3d::Zero();

	ins >> pos.x() >> pos.y() >> pos.z();
	//pos[0] = pos[0]/2.0;
	//pos[1] = pos[1]/2.0 + 2;
	//pos[2] = pos[2]/2.0;
	bbMin[0] = std::min(bbMin[0], pos.x());
	bbMin[1] = std::min(bbMin[1], pos.y());
	bbMin[2] = std::min(bbMin[2], pos.z());
	bbMax[0] = std::max(bbMax[0], pos.x());
	bbMax[1] = std::max(bbMax[1], pos.y());
	bbMax[2] = std::max(bbMax[2], pos.z());
  }

  std::cout << "total particles now: " << particles.size() << std::endl;
	std::cout << "bounding box: [" << bbMin[0] << ", " << bbMin[1] << ", "<< bbMin[2];
	std::cout << "] x [" << bbMax[0] << ", " << bbMax[1] << ", "<< bbMax[2] << "]" << std::endl;;
}

void World::saveParticleFile(const std::string& _filename) const{
  std::ofstream outs(_filename);
  
  Eigen::Vector3d pos;
	for (auto& p : particles) {
		outs << p.position.x() << " " 
			 << p.position.y() << " "
			 << p.position.z() << std::endl;
	}
	outs.close();
}

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
	  auto worldCOM = computeNeighborhoodCOM(cluster);
	  Eigen::Vector3d clusterVelocity = 
		
		computeClusterVelocity(cluster);
	  
	  Eigen::Matrix3d init;
	  init.setZero();
		

	  Eigen::Matrix3d Apq = computeApq(cluster, init, worldCOM);
	  Eigen::Matrix3d A = Apq*cluster.aInv;
	  if (nu > 0.0) A = A*cluster.Fp.inverse(); // plasticity
	  
	  //do the SVD here so we can handle fracture stuff
	  Eigen::JacobiSVD<Eigen::Matrix3d> solver(A, 
		  Eigen::ComputeFullU | Eigen::ComputeFullV);
	  
	  Eigen::Matrix3d U = solver.matrixU(), V = solver.matrixV();
	  Eigen::Vector3d sigma = solver.singularValues();
	  
	  //std::cout << "sigma " << sigma << std::endl;

	  if(fabs(sigma(0) - 1.0) > cluster.toughness){
		//if(cluster.renderWidth > toughness*cluster.width){ //doesn't improve anything
		potentialSplits.emplace_back(en.first, 
			//cluster.renderWidth - toughness*cluster.width, 
			sigma(0) - cluster.toughness, 
			V.col(0));
		//eigenvecs of S part of RS is V
	  }


	  Eigen::Matrix3d T = U*V.transpose();
	  if (nu > 0.0) T = T*cluster.Fp; // plasticity
	  
	  //auto pr = utils::polarDecomp(A);
	  
	  for(int i=0; i<cluster.neighbors.size(); i++){
		int &n = cluster.neighbors[i];
		particles[n].goalPosition += 
		  cluster.weights[i]*(T*(particles[n].restPosition - cluster.restCom) + worldCOM);
		particles[n].goalVelocity += cluster.weights[i]*clusterVelocity;
	  }
	  
	  {
	  auto timer = prof.timeName("plasticity");
	  // plasticity
	  cluster.FpNew = cluster.Fp;
	  if (nu > 0.0) {
		if (sigma(2) >= 1e-4) { // adam says: the second clause is a quick hack to avoid plasticity when sigma is degenerate
		  Eigen::Vector3d FpHat = sigma;
		  //std::cout<<FpHat(0)<<" "<<FpHat(1)<<" "<<FpHat(2)<<" => ";
		  FpHat *= 1.0/cbrt(FpHat(0) * FpHat(1) * FpHat(2));
		  //std::cout<<FpHat(0)<<" "<<FpHat(1)<<" "<<FpHat(2)<<std::endl;
		  double norm = sqrt(sqr(FpHat(0)-1.0) + sqr(FpHat(1)-1.0) + sqr(FpHat(2)-1.0));
		  double local_yield = yield + hardening * cluster.cstrain;
		  if (norm > local_yield) {	
			double gamma = std::min(1.0, nu * (norm - local_yield) / norm);
			FpHat(0) = pow(FpHat(0), gamma);
			FpHat(1) = pow(FpHat(1), gamma);
			FpHat(2) = pow(FpHat(2), gamma);
			// update cluster.Fp
			cluster.FpNew = FpHat.asDiagonal() * V.transpose() * cluster.Fp * V.determinant();
		  } 
		}
		cluster.cstrain += sqrt(sqr(sigma(0)-1.0) + sqr(sigma(1)-1.0) + sqr(sigma(2)-1.0));
	  }
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
	
	assertFinite();

  }
  
  doFracture(std::move(potentialSplits));
   
  /*
  for(auto&& en : benlib::enumerate(clusters)){
	auto& c = en.second;
	Eigen::Vector3d worldCOM = computeNeighborhoodCOM(c);
	bool updateCluster = false;
	for(auto i=0; i<c.neighbors.size(); i++) {
	  auto n = c.neighbors[i];
	  auto &p = particles[n];
	  if ((p.position - worldCOM).norm() > (1.0 + gamma) * (p.restPosition - c.restCom).norm()) {
		// create duplicate particle
		Particle q(p);
		q.clusters.clear();
		q.clusters.push_back(en.first);
		q.numClusters = 1;
		particles.push_back(q);
		c.neighbors[i] = particles.size()-1;
		// delete particle
		// remove from cluster
		// c.neighbors.erase(std::remove(c.neighbors.begin(), c.neighbors.end(), n), c.neighbors.end());
		//remove cluster from this particle
		p.clusters.erase(std::remove(p.clusters.begin(), p.clusters.end(), en.first), p.clusters.end());
		updateCluster = true;
		std::cout<<"removed an outlier "<<(particles[n].position - worldCOM).norm()<<" > "<< (1.0+gamma) * (particles[n].restPosition - c.restCom).norm()<<std::endl;
	  }
	} 
	// could update the cluster, but see below...
	}*/
  

  //cull small clusters
  auto sizeBefore = clusters.size();
  clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
		  [](const Cluster& c){
			return c.neighbors.size() < 4;
		  }), 
	  clusters.end());
  if(clusters.size() != sizeBefore){
	std::cout << "deleted " << sizeBefore - clusters.size() << " clusters" << std::endl;
  }
  updateClusterProperties(range(clusters.size()));



  
  bounceOutOfPlanes();
  for(auto&& en : benlib::enumerate(clusters)){
	auto& cluster = en.second;
	cluster.Fp = cluster.FpNew;
  }
  elapsedTime += dt;
  //std::cout << "elapsed time: " << elapsedTime << std::endl;
  for(auto& c : clusters){
	c.renderWidth = 0;
	c.worldCom = computeNeighborhoodCOM(c);
	for(auto& n : c.neighbors){
	  c.renderWidth = std::max(c.renderWidth, 
		  (c.worldCom - particles[n].position).norm());
	}
  }
  //selfCollisions();

  //printCOM();
  //std::cout<<elapsedTime<<std::endl;
}

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

  //do normal plane bounces on the backside of each plane.
  //JAL commented this out, it causes weird offsetting
//  for(auto& movingPlane : movingPlanes){
//	for(auto& p : particles){
//      //if not being pushed along outside of any plane, check for a
//      //normal bounce off the backside of the plane
//      if (!p.outsideSomeMovingPlane) {
//         movingPlane.backsideReflectBounceParticle(p, elapsedTime, epsilon);
//      }
//	}
//  }

   //handle twisting planes
  for(auto& twistingPlane : twistingPlanes){
	for(auto& p : particles){
   	  twistingPlane.twistParticle(p, elapsedTime);
   }
  }

  //do normal plane bounces on the backside of each plane.
  //JAL commented this out because the twisting planes are finite...the code is
  //there but it's slowish...
//  for(auto& twistingPlane : twistingPlanes){
//	for(auto& p : particles){
//      //if not being pushed along outside of any plane, check for a
//      //normal bounce off the backside of the plane
//      if (!p.outsideSomeMovingPlane) {
//         twistingPlane.backsideReflectBounceParticle(p, elapsedTime, epsilon);
//      }
//	}
//  }

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

void World::selfCollisions() {
  for (auto && en1 : benlib::enumerate(clusters)) {
	auto &c = en1.second;
	for (auto && en2 : benlib::enumerate(clusters)) {
	  if (en1.first == en2.first) continue;
	  auto &d = en2.second;
	  if ((c.worldCom - d.worldCom).squaredNorm() < sqr(c.width + d.width)) {
		for (auto& i : c.neighbors){
		  auto &p = particles[i];
		  if ((p.position - d.worldCom).squaredNorm() < 0.9*sqr(d.width)) {
			auto it = find(p.clusters.begin(), p.clusters.end(), en2.first);
			if (it == p.clusters.end()) {
			  //std::cout<<p.position<<std::endl<<d.worldCom<<std::endl<<d.width<<" ";
			  //std::cout<<(p.position - d.worldCom).squaredNorm() << " "<<d.width*d.width<<" "<<neighborRadius<<std::endl<<std::endl<<std::endl;;
			  //p.position = d.worldCom + d.width * (p.position - d.worldCom).normalized();
			  std::cout<<i<<" "<<en1.first<<" "<<en2.first<<" "<<(p.position - d.worldCom).norm()<<" "<<(p.restPosition - d.restCom).norm()<<std::endl;
			  //for (auto &foo : d.neighbors) {
			  //std::cout<<foo<<" ";
			  //if (foo == i) std::cout<<"maps are bad"<<std::endl;
			  //}
			  std::cout<<std::endl;
			  p.position = d.worldCom + d.width * (p.position - d.worldCom).normalized();
			}
		  }
		}
	  }
	}
  }
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

inline double cube(double x) { return x*x*x;}
//inline double poly6(double r, double h) {(r<h) ? return cube(h - r) : return 0.0;}
  

void World::makeClusters(){
  clusters.clear();
  auto r = range(particles.size());
  auto lonelyParticles = std::vector<size_t>(r.begin(), r.end());
	
  for(auto& p : particles){p.numClusters = 0;}
  
  while(!lonelyParticles.empty()) {
	auto currentParticle = lonelyParticles.back();
	lonelyParticles.pop_back();
	Cluster c;
	c.restCom = particles[currentParticle].restPosition;
	c.neighbors = restPositionGrid.getNearestNeighbors(particles,
													   c.restCom,
													   neighborRadius);
	for(auto n : c.neighbors){ ++(particles[n].numClusters);}
	clusters.push_back(c);
		
	auto it = std::remove_if(lonelyParticles.begin(), lonelyParticles.end(),
							 [&c](size_t n){
							   //neighbors aren't lonely anymore
							   return std::find(c.neighbors.begin(),
												c.neighbors.end(),
												n) != c.neighbors.end();
							 });
	lonelyParticles.erase(it, lonelyParticles.end());
  } 
	
  if (nClusters > 1) {
	std::random_device rd;
	std::mt19937 g(rd());
	// initialize cluster centers to random paricles
	while (clusters.size() < nClusters) {
	  Cluster c;
	  c.restCom = particles[g() % particles.size()].restPosition;
	  clusters.push_back(c);
	}

	/*
	{
	  for (auto& c : clusters) {
		c.neighbors.clear();
		c.weights.clear();
	  }
	  
	  for (auto i=0; i<particles.size(); i++) {
		auto& p = particles[i];
		p.clusters.clear();
		int bestNorm = std::numeric_limits<double>::infinity();
		for (auto j = 0; j<clusters.size(); j++) {
		  double newNorm = (clusters[j].restCom - p.restPosition).squaredNorm();
		  if (newNorm < bestNorm) {
			bestCluster = j;
			bestNorm = newNorm;
		  }
		}
		clusters[bestCluster].neighbors.push_back(i);
		clusters[bestCluster].weights.push_back(1.0);
		p.numClusters = 1;
		p.clusters.push_back(bestCluster);
	  }
	}
	*/	
	// fuzzy c-means loop
  	bool converged = false;
	int iters = 0;
	double sqrNeighborRadius = neighborRadius*neighborRadius;
	while (!converged || iters < 5) {
	  std::cout<<iters<<std::endl;
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
		c.neighbors = restPositionGrid.getNearestNeighbors(particles, c.restCom, neighborRadius);
		if (oldNeighbors != c.neighbors.size()) converged = false;
		c.weights.resize(c.neighbors.size());
		for (auto i=0; i<c.neighbors.size(); i++) {
		  auto n = c.neighbors[i];
		  Particle &p = particles[n];
		  double norm = (c.restCom - p.restPosition).squaredNorm();
		  double w = cube(sqrNeighborRadius-norm);
		  c.weights[i] = w;
		  p.totalweight += w;
		  p.clusters.push_back(j);
		  p.numClusters++;
		}
	  }

	  for (auto& c : clusters) {
		c.mass = 0.0;
		c.restCom = Eigen::Vector3d::Zero();
		for (auto i=0; i<c.neighbors.size(); i++) {
		  auto n = c.neighbors[i];
		  auto &p = particles[n];
		  c.weights[i] /= p.totalweight;
		  c.restCom += c.weights[i] * p.mass * p.position;
		  c.mass += c.weights[i] * p.mass;
		}
		c.restCom /= c.mass;
	  }
	} 
	std::cout<<"kmeans clustering converged in "<<iters<<std::endl;
  }
  
  updateClusterProperties(benlib::range(clusters.size()));
  
  for (auto& p : particles) {
	if (p.numClusters == 0) {
	  std::cout<<"Particle has no cluster"<<std::endl;
	  exit(0);
	}
  }
  for (auto& c : clusters) {
	if (c.mass < 1e-5) {
	  std::cout<<"Cluster has mass "<<c.mass<<" and position: "<<c.restCom<<std::endl;
	  exit(0);
	}
  }
  std::cout << "numClusters: " << clusters.size() << std::endl;

  for(auto& c : clusters){ 
	bool crossingPlane = false;
	for(auto& plane : movingPlanes){
	  bool firstSide = particles[c.neighbors.front()].position.dot(plane.normal) > plane.offset;
	  
	  for(auto n : c.neighbors){
		bool thisSide = particles[n].position.dot(plane.normal) > plane.offset;
		if(thisSide != firstSide){
		  crossingPlane = true;
		  break;
		}
	  }
	  if(crossingPlane){break;}
	}
   for(auto& plane : twistingPlanes){
	  if(crossingPlane){break;}
	  bool firstSide = particles[c.neighbors.front()].position.dot(plane.normal) > plane.offset;
	  
	  for(auto n : c.neighbors){
		bool thisSide = particles[n].position.dot(plane.normal) > plane.offset;
		if(thisSide != firstSide){
		  crossingPlane = true;
		  break;
		}
	  }
	}
   for(auto& plane : tiltingPlanes){
	  if(crossingPlane){break;}
	  bool firstSide = particles[c.neighbors.front()].position.dot(plane.normal) > plane.offset;
	  
	  for(auto n : c.neighbors){
		bool thisSide = particles[n].position.dot(plane.normal) > plane.offset;
		if(thisSide != firstSide){
		  crossingPlane = true;
		  break;
		}
	  }
	}
	if(crossingPlane){
	  c.toughness = 10*toughness;//std::numeric_limits<double>::infinity();
	} else {
	  c.toughness = toughness;
	}
  }

  


}


void World::strainLimitingIteration(){
  const double gammaSquared = gamma * gamma; 
  for(auto& p : particles){
		p.goalPosition.setZero();
  }
  for(auto& c : clusters){
	Eigen::Vector3d worldCOM = computeNeighborhoodCOM(c);
	Eigen::Matrix3d init; 
	init.setZero();
	
	Eigen::Matrix3d Apq = computeApq(c, init, worldCOM);
	Eigen::Matrix3d A = Apq*c.aInv;
	if (nu > 0.0) A = A*c.Fp.inverse(); // plasticity
	auto pr = utils::polarDecomp(A);

	Eigen::Matrix3d T = pr.first;
	if (nu > 0.0) T = T * c.Fp;

	for (int i=0; i<c.neighbors.size(); i++) {
	  //for(auto n : c.neighbors){
	  int &n = c.neighbors[i];
	  auto &q = particles[n];
	  Eigen::Vector3d rest = (q.restPosition - c.restCom);
	  Eigen::Vector3d goal = T*(rest) + worldCOM;
	  double ratio = (goal-q.position).squaredNorm() / 
		(c.width*c.width);
	  
	  if (ratio > gammaSquared) {
		q.goalPosition += //(1.0/p.clusterMass)*
		  c.weights[i] * (goal + 
			  sqrt(gammaSquared/ratio) * 
			  (q.position - goal));
	  } else {
		q.goalPosition += //(1.0/p.clusterMass)*
		  c.weights[i] * q.position;
	  }
	}
	
  }
  
  for(auto& p : particles){
	if(p.numClusters > 0){
	  //p.goalPosition /= p.numClusters;
	} else {
	  p.goalPosition = p.position;
	}
	p.position = omega*p.goalPosition + (1.0-omega)*p.position;
  }
}

// this will be wrong after particles are duplicated
void World::printCOM() const{
  Eigen::Vector3d worldCOM = 
	std::accumulate(particles.begin(), particles.end(),
					Eigen::Vector3d{0,0,0},
					[](Eigen::Vector3d acc, const Particle& p){
					  return acc + p.mass*p.position;
					});
  double totalMass = std::accumulate(particles.begin(), particles.end(),
									 0.0, 
									 [](double acc, const Particle& p){
									   return acc + p.mass;
									 });
  
  worldCOM /= totalMass;
  std::cout << worldCOM.x() << std::endl;;
}


Eigen::Vector3d World::computeNeighborhoodCOM(const Cluster& c) const {
  return sumWeightedWorldCOM(c.neighbors, c.weights);
  //positions weighted by (mass/numClusters)
  /*
  return std::accumulate(c.neighbors.begin(), c.neighbors.end(),
	  Eigen::Vector3d(0.0, 0.0, 0.0),
	  [this](Eigen::Vector3d acc, int n){
		return acc + (particles[n].mass/particles[n].numClusters)*
		  particles[n].position;
	  })/sumWeightedMass(c.neighbors, c.weights);
  */
}

// adam says: this should take weights into account
Eigen::Matrix3d World::computeApq(const Cluster& c, 
								  const Eigen::Matrix3d& init,
								  const Eigen::Vector3d& worldCOM) const{
  return std::accumulate(c.neighbors.begin(), c.neighbors.end(),
						 init,
						 [this,&c, &worldCOM]
						 (const Eigen::Matrix3d& acc, int n) -> 
						 Eigen::Matrix3d{
						   Eigen::Vector3d pj = particles[n].position - worldCOM;
						   Eigen::Vector3d qj = particles[n].restPosition - c.restCom;
						   //return acc + (particles[n].mass/particles[n].numClusters)*
						   return acc + (particles[n].mass)*
							 pj*qj.transpose();
						 });
} 	

Eigen::Vector3d World::computeClusterVelocity(const Cluster &c) const {//(const std::vector<int> &indices, const std::vector<double> &weights) const {
	double mass = 0.0;
	Eigen::Vector3d vel = Eigen::Vector3d::Zero();
	std::vector<int>::const_iterator ii = c.neighbors.begin();
	std::vector<double>::const_iterator di = c.weights.begin();
	for (; ii!=c.neighbors.end(); ii++, di++) {
	  mass += (*di)*particles[*ii].mass;
	  vel += (*di)*particles[*ii].mass*particles[*ii].velocity;
	}
	return (vel / mass);
  }

/*
Eigen::Vector3d World::computeClusterVelocity(const Cluster& c) const {
  return 
	std::accumulate(c.neighbors.begin(),
					c.neighbors.end(),
					Eigen::Vector3d(0.0, 0.0, 0.0),
					[this](Eigen::Vector3d acc, int n){
					  return acc + (particles[n].mass/
									particles[n].numClusters)*
						particles[n].velocity;
					})/
	std::accumulate(c.neighbors.begin(), c.neighbors.end(),
					0.0,
					[this](double acc, int n){
					  return acc + 
						particles[n].mass/particles[n].numClusters;
					});
  
					}*/

void World::countClusters(){
  for(auto& p : particles){
	p.numClusters = 0;
	p.clusters.clear();
  }
  for(auto cInd : benlib::range(clusters.size())){
	for(auto& i : clusters[cInd].neighbors){
	  ++(particles[i].numClusters);
	  particles[i].clusters.push_back(cInd);
	}
  }
  for(auto& p : particles){assert(p.numClusters == p.clusters.size());}
}

template <typename Container>
void World::updateClusterProperties(const Container& clusterIndices){
  countClusters(); //TODO, just touch a subset...

  // compute cluster mass, com, width, and aInv
  for(auto cIndex : clusterIndices){
	auto& c = clusters[cIndex];
	//c.Fp.setIdentity(); // plasticity
	c.mass = sumWeightedMass(c.neighbors, c.weights);
	assert(c.mass >= 0);
	c.restCom = sumWeightedRestCOM(c.neighbors, c.weights);
	c.worldCom = computeNeighborhoodCOM(c);

	if(!c.restCom.allFinite()){std::cout << c.restCom << std::endl;}
	if(!c.worldCom.allFinite()){std::cout << c.worldCom << std::endl;}

	assert(c.restCom.allFinite());
	assert(c.worldCom.allFinite());
	c.width = 0.0;
	for(auto n : c.neighbors){
	  c.width = std::max(c.width, (c.restCom - particles[n].restPosition).norm());
	} 
	c.renderWidth = c.width; //it'll get updated soon enough
	//assert(c.width >= 0);
	
	c.aInv.setZero();  
	// adam says: this should take weights into account
	c.aInv = 
	  std::accumulate(c.neighbors.begin(), c.neighbors.end(),
		  c.aInv,
		  [this, &c]
		  (const Eigen::Matrix3d& acc, int n) -> 
		  Eigen::Matrix3d {
			Eigen::Vector3d qj = particles[n].restPosition - c.restCom;
			return acc + (particles[n].mass)*
			  qj*qj.transpose();
		  });
	
	//do pseudoinverse
	Eigen::JacobiSVD<Eigen::Matrix3d> solver(c.aInv, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Vector3d sigInv;
	for(auto i : range(3)){
	  if (solver.singularValues()(i) < 1e-6) c.Fp.setIdentity();
	  sigInv(i) = fabs(solver.singularValues()(i)) > 1e-6 ? 1.0/solver.singularValues()(i) : 0;
	}
	c.aInv = solver.matrixV()*sigInv.asDiagonal()*solver.matrixU().transpose();//c.aInv.inverse().eval();
	if(!c.aInv.allFinite()){
	  std::cout << c.aInv << std::endl;
	  std::cout << solver.singularValues() << std::endl;
	}
	assert(c.aInv.allFinite());
  }

}


void World::doFracture(std::vector<World::FractureInfo> potentialSplits){
  auto timer = prof.timeName("fracture");
  //do fracture
  std::sort(potentialSplits.begin(), potentialSplits.end(),
	  [](const FractureInfo& a, const FractureInfo& b){
		return std::get<1>(a) < std::get<1>(b);
	  });
  int count = 0;
  if(!potentialSplits.empty()){
	std::cout << "potential splits: " << potentialSplits.size() << std::endl;
  }
  bool aCancel = false;
  for(auto &ps : potentialSplits){
	//if (++count > 10) break;
	size_t cIndex = std::get<0>(ps);

	auto worldCOM = clusters[cIndex].worldCom;

	//recompute A matrix:
	Eigen::Matrix3d init;
	init.setZero();

	Eigen::Matrix3d Apq = computeApq(clusters[cIndex], init, worldCOM);
	Eigen::Matrix3d A = Apq*clusters[cIndex].aInv;
	if (nu > 0.0) A = A*clusters[cIndex].Fp.inverse(); // plasticity
	
	//do the SVD here so we can handle fracture stuff
	Eigen::JacobiSVD<Eigen::Matrix3d> solver(A, Eigen::ComputeFullV);
	
	Eigen::Vector3d sigma = solver.singularValues();
	if(fabs(sigma(0) - 1) < clusters[cIndex].toughness){
	  if(!aCancel){
		aCancel = true;
		std::cout << "cancelled fracture with updated stuff" << std::endl;
	  }
	  continue;
	}
	
	Eigen::Vector3d splitDirection = solver.matrixV().col(0);

	//doesn't work... 
	//just erase the cluster
	//clusters.erase(clusters.begin() + cIndex);
	//updateClusterProperties();
	//break;


	

	//auto& cluster = clusters[cIndex];	  
	// adam says: why is this a bad idea?  clusters[cIndex] is ugly and shows up a lot.
	//ben says: when I push_back, the reference gets invalidated if the vector reallocates (which bit me).


	//if(cluster.neighbors.size() < 10){ continue;}

	auto it = std::partition(clusters[cIndex].neighbors.begin(),
		clusters[cIndex].neighbors.end(),
		[&worldCOM, &splitDirection, this](int ind){
		  //which side of the split is it on?
		  return (worldCOM - particles[ind].position).dot(splitDirection) > 0;
		});
	auto oldSize = std::distance(clusters[cIndex].neighbors.begin(), it);
	auto newSize = std::distance(it, clusters[cIndex].neighbors.end());
	if(newSize == 0 || oldSize == 0){ continue;}
	//if(oldSize < 4 || newSize < 4){ continue;}

	//expected to be mostly in teh x direction for the cube example, and it was
	//std::cout << "split direction: " << splitDirection << std::endl;
	
	//make a new cluster
	Cluster newCluster;
	newCluster.neighbors.assign(it, clusters[cIndex].neighbors.end());

	// copy relevant variables
	newCluster.Fp = clusters[cIndex].Fp; // plasticity
	newCluster.FpNew = clusters[cIndex].FpNew;
	newCluster.cstrain = clusters[cIndex].cstrain; // plasticity
	newCluster.toughness = clusters[cIndex].toughness;
	// we will want to copy toughness here as well...
	
	//delete the particles from the old one
	clusters[cIndex].neighbors.erase(clusters[cIndex].neighbors.begin() + oldSize, 
		clusters[cIndex].neighbors.end());
	
	clusters.push_back(newCluster);	  
	
	updateClusterProperties(std::initializer_list<size_t>{cIndex, clusters.size() -1});
	//std::cout << "numClusters: " << clusters.size() << std::endl;
	
	//std::cout << "min cluster size: " << std::min_element(clusters.begin(), clusters.end(),
	//		[](const Cluster& a, const Cluster& b){
	//		  return a.neighbors.size() < b.neighbors.size();})->neighbors.size() << std::endl;
	
	
	{
	  auto timerTwo = prof.timeName("propagate");
	  //split from other clusters
	  std::vector<int> allParticles(clusters[cIndex].neighbors.size() + newCluster.neighbors.size());
	  std::copy(newCluster.neighbors.begin(), newCluster.neighbors.end(),
		  std::copy(clusters[cIndex].neighbors.begin(), clusters[cIndex].neighbors.end(), allParticles.begin()));
	  std::vector<size_t> affectedClusters; //keep sorted
	  for(auto& member : allParticles){
		auto& particle = particles[member];
		for(auto thisIndex : particle.clusters){
		  //insert into sorted array
		  auto it = std::lower_bound(affectedClusters.begin(), affectedClusters.end(), thisIndex);
		  if(it == affectedClusters.end() || *it != thisIndex){
			affectedClusters.insert(it, thisIndex);
		  }
		  
		  auto& thisCluster = clusters[thisIndex];
		  //auto thisClusterCOM = computeNeighborhoodCOM(thisCluster);
		  if(((particle.position - worldCOM).dot(splitDirection) >= 0) !=
			  ((thisCluster.worldCom - worldCOM).dot(splitDirection) >= 0 )){
			//remove from cluster
			thisCluster.neighbors.erase(
				std::remove(thisCluster.neighbors.begin(),
					//thisCluster.neighbors.end(), thisIndex), thisCluster.neighbors.end());
					thisCluster.neighbors.end(), member), thisCluster.neighbors.end());
			//remove cluster from this
			particle.clusters.erase(
				std::remove(particle.clusters.begin(), particle.clusters.end(),
					thisIndex), particle.clusters.end());
			
		  }
		}
	  }
	  updateClusterProperties(affectedClusters);
	  //	for(auto& c : affectedClusters){
	  //	  clusters[c].toughness *= 0.995;
	//	}
	
	//break;
	}
  }
}	

void World::dumpParticlePositions(const std::string& filename) const{
  std::ofstream outs(filename, std::ios_base::binary | std::ios_base::out);
  size_t numParticles = particles.size();
  outs.write(reinterpret_cast<const char*>(&numParticles), sizeof(numParticles));
  std::vector<float> positions(3*numParticles);
  for(auto i : range(particles.size())){
	positions[3*i    ] = particles[i].position.x();
	positions[3*i + 1] = particles[i].position.y();
	positions[3*i + 2] = particles[i].position.z();
  }
  outs.write(reinterpret_cast<const char*>(positions.data()), 
	  3*numParticles*sizeof(typename decltype(positions)::value_type));

}

void World::dumpColors(const std::string& filename) const {
  std::ofstream outs(filename, std::ios_base::binary | std::ios_base::out);
  size_t numParticles = particles.size();
  outs.write(reinterpret_cast<const char*>(&numParticles), sizeof(numParticles));
  std::vector<float> colors(3*numParticles);
  for(auto&& pr : benlib::enumerate(particles)) {
	//find nearest cluster
	auto i = pr.first;
	auto& p = pr.second;

	int min_cluster = p.clusters[0];
	auto com = clusters[min_cluster].worldCom;
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
	if(clusters[min_cluster].neighbors.size() > 1) {
	  //      sqrt(min_sqdist) < 0.55*clusters[min_cluster].renderWidth) {
	  //glColor4d(rgb.r, rgb.g, rgb.b, 0.8);
	  colors[3*i]  = rgb.r;
	  colors[3*i+ 1]  = rgb.g;
	  colors[3*i+ 2]  = rgb.b;
	} else {
	  colors[3*i] = 1;
	  colors[3*i + 1] = 1;
	  colors[3*i + 2] = 1;
	}
  }
  outs.write(reinterpret_cast<const char*>(colors.data()),
	  3*numParticles*sizeof(float));
}
