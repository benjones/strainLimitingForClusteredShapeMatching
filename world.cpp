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

#include <random>

#include "enumerate.hpp"
using benlib::enumerate;
#include "range.hpp"
using benlib::range;

void World::draw(SDL_Window* window) const {

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



  //draw planes
  glBegin(GL_QUADS);
  for(auto&& pr : enumerate(planes)){
	const auto i = pr.first;
	const auto& plane = pr.second;
	Eigen::Vector3d tangent1, tangent2;

	const Eigen::Vector3d norm = plane.head(3);
	const double w = plane.w();

	tangent1 = norm.cross(Eigen::Vector3d{1,0,0});
	if(tangent1.isZero(1e-3)){
	  tangent1 = norm.cross(Eigen::Vector3d{0,0,1});
	  if(tangent1.isZero(1e-3)){
		tangent1 = norm.cross(Eigen::Vector3d{0,1,0});
	  }
	}
	tangent1.normalize();

	tangent2 = norm.cross(tangent1);
	tangent2.normalize(); //probably not necessary
	
	const double sos = norm.dot(norm);
	const Eigen::Vector3d supportPoint{norm.x()*w/sos,
		norm.y()*w/sos,
		norm.z()*w/sos};


	glColor4d(0.5, static_cast<double>(i)/planes.size(),
			  0.5, 1);

	const double size = 100;
	glNormal3dv(norm.data());
	glVertex3dv((supportPoint + size*(tangent1 + tangent2)).eval().data());
	glVertex3dv((supportPoint + size*(-tangent1 + tangent2)).eval().data());
	glVertex3dv((supportPoint + size*(-tangent1 - tangent2)).eval().data());
	glVertex3dv((supportPoint + size*(tangent1  - tangent2)).eval().data());
  }
  glEnd();


  
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



void World::loadFromJson(const std::string& _filename){

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

  auto gravityIn = root["gravity"];
  if(!gravityIn.isNull() && gravityIn.isArray() && gravityIn.size() == 3){
	gravity.x() = gravityIn[0].asDouble();
	gravity.y() = gravityIn[1].asDouble();
	gravity.z() = gravityIn[2].asDouble();
  } else {
	gravity = Eigen::Vector3d{0, -9.81, 0};
  }

  auto planesIn =  root["planes"];
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
  
  double mass = root.get("mass", 0.1).asDouble();
  for(auto& p : particles){ p.mass = mass;}

  
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

void World::timestep(){
  
  //scope block for profiler
  {
	auto timer = prof.timeName("shape matching");
	for(auto& p : particles){
	  p.oldPosition = p.position;
	  p.goalPosition.setZero();
	  p.goalVelocity.setZero();
	}
	
	for(auto& c : clusters){
	  auto worldCOM = computeNeighborhoodCOM(c);
	  Eigen::Vector3d clusterVelocity = 

		computeClusterVelocity(c);

		Eigen::Matrix3d init;
		init.setZero();
		
		Eigen::Matrix3d Apq = computeApq(c, init, worldCOM);
		Eigen::Matrix3d A = Apq*c.aInv;
		auto pr = utils::polarDecomp(A);
		
		for(auto n : c.neighbors){
		  particles[n].goalPosition += 
			(pr.first*(particles[n].restPosition - c.restCom) + worldCOM);
		  particles[n].goalVelocity += clusterVelocity;
		}
	}
	
	for(auto& p : particles){
	  p.goalPosition /= p.numClusters;
	  p.goalVelocity /= p.numClusters;
	  p.velocity += dt * gravity + (alpha/dt)*(p.goalPosition- p.position) + 
		springDamping*(p.goalVelocity - p.velocity); 
	  p.position += dt * p.velocity;
	}
	
  }

  //scope block for profiler
  {
	auto timer = prof.timeName("strainLimiting");
	for(auto iter : range(numConstraintIters)){
	  (void)iter; //unused
	  strainLimitingIteration();
	}
	for(auto& p : particles){
	  p.velocity = (1.0/dt)*(p.position - p.oldPosition);
	}
	
  }
  
  bounceOutOfPlanes();
  
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
	// we'll use numClusters to keep track of the best cluster
	for (auto& p : particles) p.numClusters = -1;
		
	// kmeans loop
	bool converged = false;
	int iters = 0;
	while (!converged) {
	  converged = true;
	  iters++;
	  for (auto& c : clusters) c.neighbors.clear();
	  
	  for (auto i=0; i<particles.size(); i++) {
		auto& p = particles[i];
		int bestCluster = 0;
		double bestNorm = (clusters[0].restCom - p.restPosition).squaredNorm();
		for (auto j = 0; j<clusters.size(); j++) {
		  double newNorm = (clusters[j].restCom - p.restPosition).squaredNorm();
		  if (newNorm < bestNorm) {
			bestCluster = j;
			bestNorm = newNorm;
		  }
		}
		clusters[bestCluster].neighbors.push_back(i);
		if (p.numClusters != bestCluster) converged = false;
		p.numClusters = bestCluster;
	  }
	  
	  for (auto& c : clusters) {
		double mass = std::accumulate(c.neighbors.begin(), c.neighbors.end(),
									  0.0,
									  [this](double acc, int n){
										return acc + 
										particles[n].mass;
									  });
		
		if (mass > 1e-5) 
		  c.restCom = std::accumulate(c.neighbors.begin(), c.neighbors.end(),
									  Eigen::Vector3d(0.0, 0.0, 0.0),
									  [this](Eigen::Vector3d acc, int n){
										return acc + particles[n].mass*
										particles[n].position;
									  }) / mass;
	  }
	}	
	
	// count numClusters and initialize neighborhoods
	for (auto& p : particles) p.numClusters = 0;
	for (auto& c : clusters) {
	  c.neighbors = restPositionGrid.getNearestNeighbors(particles, c.restCom, neighborRadius);
	  for(auto n : c.neighbors) ++(particles[n].numClusters);
	}
	std::cout<<"kmeans clustering converged in "<<iters<<std::endl;
  }
  
  // compute cluster mass, com, width, and aInv
  for (auto& c : clusters) {
	c.mass = 	std::accumulate(c.neighbors.begin(),
								c.neighbors.end(),
								0.0,
								[this](double acc, int n){
								  return acc + particles[n].mass / particles[n].numClusters;
								});
	c.restCom = computeNeighborhoodCOM(c);
	c.width = 0.0;
	for(auto n : c.neighbors){
	  c.width = std::max(c.width, (c.restCom - particles[n].restPosition).norm());
	} 

	c.aInv.setZero();  
	c.aInv = 
	  std::accumulate(c.neighbors.begin(), c.neighbors.end(),
					  c.aInv,
					  [this, &c]
					  (const Eigen::Matrix3d& acc, int n) -> 
					  Eigen::Matrix3d {
						Eigen::Vector3d qj = particles[n].restPosition - c.restCom;
						// return acc + (particles[n].mass/particles[n].numClusters)*
						return acc + (particles[n].mass)*
						  qj*qj.transpose();
					  });
	
	c.aInv = c.aInv.inverse().eval();
  }

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
	auto pr = utils::polarDecomp(A);
	
	for(auto n : c.neighbors){
	  auto &q = particles[n];
	  Eigen::Vector3d rest = (q.restPosition - c.restCom);
	  Eigen::Vector3d goal = pr.first*(rest) + worldCOM;
	  double ratio = (goal-q.position).squaredNorm() / 
		(c.width*c.width);
	  
	  if (ratio > gammaSquared) {
		q.goalPosition += //(1.0/p.clusterMass)*
		  (goal + 
		   sqrt(gammaSquared/ratio) * 
		   (q.position - goal));
	  } else {
		q.goalPosition += //(1.0/p.clusterMass)*
		  q.position;
	  }
	}
  }
  
  for(auto& p : particles){
		p.goalPosition /= p.numClusters;
		p.position = omega*p.goalPosition + (1.0-omega)*p.position;
  }
}

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
  //positions weighted by (mass/numClusters)
  return std::accumulate(c.neighbors.begin(), c.neighbors.end(),
						 Eigen::Vector3d(0.0, 0.0, 0.0),
						 [this](Eigen::Vector3d acc, int n){
						   return acc + (particles[n].mass/particles[n].numClusters)*
							 particles[n].position;
						 })/
	std::accumulate(c.neighbors.begin(), c.neighbors.end(),
					0.0,
					[this](double acc, int n){
					  return acc + 
						particles[n].mass/particles[n].numClusters;
					});
}

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
  
}
