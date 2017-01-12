#include "world.h"
#include <fstream>

#include "json/json.h"
#include "utils.h"
#include "color_spaces.h"

#include "enumerate.hpp"
using benlib::enumerate;
#include "range.hpp"
using benlib::range;




void World::loadFromJson(const std::string& _filename){
  elapsedTime = 0;
  filename = _filename;

  clusters.clear();
  particles.clear();
  
  std::ifstream ins(filename);
  
  Json::Value root;
  Json::Reader jReader;

  if(!jReader.parse(ins, root)){
	std::cout << "couldn't read input file: " << filename << '\n'
			  << jReader.getFormattedErrorMessages() << std::endl;
	exit(1);
  }

  clusteringParams.neighborRadius = root.get("neighborRadius", 0.1).asDouble();
  clusteringParams.nClusters = root.get("nClusters", 1).asInt();
  clusteringParams.neighborRadiusMax =
	root.get("neighborRadiusMax", std::numeric_limits<double>::max()).asDouble();
  clusteringParams.nClustersMax = root.get("nClustersMax", std::numeric_limits<int>::max()).asInt();
  clusteringParams.clusterItersMax = root.get("clusterItersMax", 10000).asInt();
  clusteringParams.clusteringAlgorithm = root.get("clusteringAlgorithm", 0).asInt();
  clusteringParams.clusterOverlap = root.get("clusterOverlap", 0.0).asDouble();
  clusteringParams.clusterKernel = root.get("clusterKernel", 0).asInt();
  clusteringParams.kernelWeight = root.get("kernelWeight", 2.0).asDouble();
  clusteringParams.blackhole = root.get("blackhole", 1.0).asDouble();

  double mass = root.get("mass", 0.1).asDouble();


  
  auto particleFilesIn = root["particleFiles"];
  for(auto i : range(particleFilesIn.size())){
	auto particleInfo = readParticleFile(particleFilesIn[i].asString());

	std::vector<Particle> newParticles;
	newParticles.reserve(particleInfo.positions.size());
	for(const auto & pos : particleInfo.positions){
	  newParticles.emplace_back();
	  newParticles.back().position = pos;
	  newParticles.back().restPosition = pos;
	  newParticles.back().velocity = Eigen::Vector3d::Zero();
	  newParticles.back().mass = mass;
	}
	auto newClusters = iterateMakeClusters(newParticles, clusteringParams);
	mergeClusters(newParticles, newClusters);

	
  }

  auto movingParticleFilesIn = root["movingParticleFiles"];
  for(auto i : range(movingParticleFilesIn.size())){
	auto& movingP = movingParticleFilesIn[i];
	auto particleInfo = readParticleFile(movingP["filename"].asString());
	Eigen::Vector3d offset = Eigen::Vector3d::Zero();
	double scale = movingP.get("scale", 1.0).asDouble();
	auto& os = movingP["offset"];
	if(!os.isNull()){
	  offset.x() = os[0].asDouble();
	  offset.y() = os[1].asDouble();
	  offset.z() = os[2].asDouble();
	}
	
	Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
	auto& vel = movingP["velocity"];
	if(!vel.isNull()){
	  velocity.x() = vel[0].asDouble();
	  velocity.y() = vel[1].asDouble();
	  velocity.z() = vel[2].asDouble();
	}
	std::vector<Particle> newParticles;
	newParticles.reserve(particleInfo.positions.size());
	for(const auto& pos : particleInfo.positions){
	  newParticles.emplace_back();
	  newParticles.back().position = scale*pos + offset;
	  newParticles.back().restPosition = scale*pos + offset;
	  newParticles.back().velocity = velocity;
	  newParticles.back().mass = mass;
	}
	
	auto newClusters = iterateMakeClusters(newParticles, clusteringParams);
	mergeClusters(newParticles, newClusters);
  }


  /* unused in all our examples, and doesn't really belong, refactored to the vis stuff --Ben
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
	}*/

  dt = root.get("dt",1/60.0).asDouble();
  numConstraintIters = root.get("numConstraintIters", 5).asInt();
  alpha = root.get("alpha", 1.0).asDouble();
  omega = root.get("omega", 1.0).asDouble();
  gamma = root.get("maxStretch", 0.1).asDouble();
  gamma = root.get("gamma", gamma).asDouble();
  springDamping = root.get("springDamping", 0.0).asDouble();
  toughness = root.get("toughness", std::numeric_limits<double>::infinity()).asDouble();
  toughnessBoost = root.get("toughnessBoost", 0.0).asDouble();
  toughnessFalloff = root.get("toughnessFalloff", std::numeric_limits<double>::infinity()).asDouble();

  // plasticity parameters
  yield = root.get("yield", 0.0).asDouble();
  nu = root.get("nu", 0.0).asDouble();
  hardening = root.get("hardening", 0.0).asDouble();
  
  collisionRestitution = root.get("collisionRestitution", 0.5).asDouble();
  collisionGeometryThreshold = root.get("collisionGeometryThreshold", 0.5).asDouble();
  std::cout<<"collisionGeometryThreshold "<<collisionGeometryThreshold<<std::endl;
  outlierThreshold = root.get("outlierThreshold", 2.0).asDouble();

  auto fractureIn = root["fracture"];
  if (!fractureIn.isNull() && fractureIn.asString().compare("off") == 0) {
	fractureOn = false;
	std::cout<<"Fracture off"<<std::endl;
  } else fractureOn = true;

  auto delayRepeatedFractureIn = root["delayRepeatedFracture"];
  if (!delayRepeatedFractureIn.isNull() && delayRepeatedFractureIn.asString().compare("off") == 0) {
	delayRepeatedFracture = false;
	std::cout<<"delayRepeatedFracture off"<<std::endl;
  } else delayRepeatedFracture = true;

  auto selfCollisionsIn = root["selfCollisions"];
  if (!selfCollisionsIn.isNull() && selfCollisionsIn.asString().compare("off") == 0) {
	selfCollisionsOn = false;
	std::cout<<"Self Collisions Off"<<std::endl;
  } else selfCollisionsOn = true;

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



  //do this in mergeClusters instead
  //for(int i=0; i<(int)particles.size(); i++) particles[i].id = i;

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


  Eigen::Vector3d initialStretch;
  auto& initialStretchIn = root["initialStretch"];
  if(!initialStretchIn.isNull() && initialStretchIn.isArray() && initialStretchIn.size() == 3){
	initialStretch.x() = initialStretchIn[0].asDouble();
	initialStretch.y() = initialStretchIn[1].asDouble();
	initialStretch.z() = initialStretchIn[2].asDouble();
  } else {
	initialStretch = Eigen::Vector3d{1.0, 1.0, 1.0};
  }


  
  


  for (auto& c : clusters) {
	Eigen::Matrix3d Aqq;
	Aqq.setZero();
	for (auto &member : c.members) {
	  auto &p = particles[member.index];
	  Eigen::Vector3d qj = p.restPosition - c.restCom;
	  Aqq += qj * qj.transpose();
	}
	Eigen::JacobiSVD<Eigen::Matrix3d> solver(Aqq, Eigen::ComputeFullU | Eigen::ComputeFullV);

	double min, max;
	Eigen::Vector3d n;
	for (unsigned int i=0; i<3; i++) {
	  n = solver.matrixV().col(i);
	  min = std::numeric_limits<double>::max();
	  max = -std::numeric_limits<double>::max();
	  for (auto &member : c.members) {
		auto &p = particles[member.index];
		double d = n.dot(p.restPosition);
		if (d < min) min = d;
		if (d > max) max = d;
	  }
	  double center = n.dot(c.cg.c);
	  //std::cout<<max - min<<" "<<neighborRadius<<std::endl;
	  //if (max - min < 1.5*neighborRadius) {
	  if (max - center < collisionGeometryThreshold*neighborRadius) {
		std::cout<<" adding plane: "<<n(0)<<"x + "<<n(1)<<"y + "<<n(2)<<"z + "<<-max<<std::endl;
		c.cg.addPlane(n, -max);
	  }
	  if (center - min < collisionGeometryThreshold*neighborRadius) {
		std::cout<<" adding plane: "<<-n(0)<<"x + "<<-n(1)<<"y + "<<-n(2)<<"z + "<<min<<std::endl;
		c.cg.addPlane(-n, min);
	  }
	}
  }

  updateClusterProperties(benlib::range(clusters.size()));

  for (auto& c : clusters) {
	if (c.mass < 1e-5) {
	  std::cout<<"Cluster has mass "<<c.mass<<" and position: "<<c.restCom<<std::endl;
	  exit(0);
	}
  }


  
  setupPlaneConstraints();

  //color the particles... hopefully the clusters have the COMS correct here
  for(auto& p : particles){
	auto closestCluster = std::min_element(
		p.clusters.begin(), p.clusters.end(),
		[&p,this](int a, int b){
		  return (p.position - clusters[a].worldCom).squaredNorm() <
		  (p.position - clusters[b].worldCom).squaredNorm();
		});
	
	p.color = HSLColor(2.0*acos(-1)*(*closestCluster%12)/12.0, 0.7, 0.7).to_rgb();

	
  }
  
  
  for (auto &p : particles) {
	p.position.x() *= initialStretch[0];
	p.position.y() *= initialStretch[1];
	p.position.z() *= initialStretch[2];
  }
  //apply initial rotation/scaling, etc
  //todo, read this from json
  //Eigen::Vector3d axis{0,0,1};
  //for (auto &p : particles) {
	//p.position.x() *= 2.5;
	//p.velocity = 0.5*p.position.cross(axis);
	//auto pos = p.position;
	//p.position[0] = 0.36*pos[0] + 0.48*pos[1] - 0.8*pos[2];
	//p.position[1] = -0.8*pos[0] + 0.6*pos[1];// - 0.8*pos[2];
	//p.position[2] = 0.48*pos[0] + 0.64*pos[1] + 0.6*pos[2];
  //}
}




World::ParticleSet World::readParticleFile(const std::string& _filename){


  std::ifstream ins(_filename);
  
  ParticleSet ret;

  Eigen::Vector3d pos;
  ins >> pos.x() >> pos.y() >> pos.z();

  ret.bbMin = pos;
  ret.bbMax = pos;
  
  while(ins){
	ret.positions.push_back(pos);
	
	ret.bbMin - ret.bbMin.cwiseMin(pos);
	ret.bbMax - ret.bbMax.cwiseMax(pos);

	ins >> pos.x() >> pos.y() >> pos.z();
  }

  std::cout << "total particles in file " << _filename << ": " << ret.positions.size() << std::endl;
  std::cout << "bounding box: [" << ret.bbMin[0] << ", " << ret.bbMin[1] << ", "<< ret.bbMin[2];
  std::cout << "] x [" << ret.bbMax[0] << ", " << ret.bbMax[1] << ", "<< ret.bbMax[2] << "]" << std::endl;
  return ret;
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
	/*	//find nearest cluster
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

	RGBColor rgb = HSLColor(2.0*acos(-1)*(min_cluster%12)/12.0, 0.7, 0.7).to_rgb();*/
	//RGBColor rgb = HSLColor(2.0*acos(-1)*min_cluster/clusters.size(), 0.7, 0.7).to_rgb();
	//if(clusters[min_cluster].members.size() > 1) {
	  //      sqrt(min_sqdist) < 0.55*clusters[min_cluster].renderWidth) {
	  //glColor4d(rgb.r, rgb.g, rgb.b, 0.8);
	auto i = pr.first;
	auto& p = pr.second;
	  colors[3*i]  = p.color.r;//rgb.r;
	  colors[3*i+ 1]  = p.color.g;//rgb.g;
	  colors[3*i+ 2]  = p.color.b;//rgb.b;
	  //always use the color...
	  /*} else {
	  colors[3*i] = 1;
	  colors[3*i + 1] = 1;
	  colors[3*i + 2] = 1;
	  }*/
  }
  outs.write(reinterpret_cast<const char*>(colors.data()),
	  3*numParticles*sizeof(float));
}


void World::dumpClippedSpheres(const std::string& filename) const {
  std::ofstream outs(filename);
  
  outs << clusters.size() << '\n';
  for(const auto& cluster : clusters){

	auto visTrans = cluster.getVisTransform();

   //Eigen::Vector3d newCenter = cluster.cg.c + visTrans.col(3).head(3);
   Eigen::Vector4d newCenter(cluster.cg.c(0), cluster.cg.c(1), cluster.cg.c(2), 1);
   newCenter = visTrans * newCenter;


	outs << cluster.cg.planes.size() << '\n'
		 << newCenter.x() << ' '
		 << newCenter.y() << ' '
		 << newCenter.z() << '\n'
		 << cluster.cg.r << '\n';
	for(const auto& plane : cluster.cg.planes){
	  //Eigen::Vector3d newNormal = visTrans.block(0,0,3,3)*plane.first;
     // outs << newNormal.x() << ' '
	//	   << newNormal.y() << ' '
	//	   << newNormal.z() << ' '
	//	   << plane.second << '\n';
	  Eigen::Vector4d norm4(plane.first(0), plane.first(1), plane.first(2), 0);
     Eigen::Vector4d pt4(cluster.cg.c(0), cluster.cg.c(1), cluster.cg.c(2), 1);
     double d = norm4.dot(pt4) + plane.second;
     pt4 = pt4 - d*norm4;

     norm4 = visTrans * norm4;
     pt4 = visTrans * pt4;

     Eigen::Vector3d norm3(norm4(0), norm4(1), norm4(2));
     Eigen::Vector3d pt3(pt4(0), pt4(1), pt4(2));
     d = norm3.dot(pt3);
      
     outs << norm3.x() << ' '
		   << norm3.y() << ' '
	      << norm3.z() << ' '
		   << -d << '\n';
      
	}


  }
}


void World::setupPlaneConstraints(){

  for(auto& c : clusters){ 
	bool crossingPlane = false;
	for(auto& plane : movingPlanes){
	  bool firstSide = particles[c.members.front().index].position.dot(plane.normal) > plane.offset;
	  
	  for(const auto& member : c.members){
		bool thisSide = particles[member.index].position.dot(plane.normal) > plane.offset;
		if(thisSide != firstSide){
		  crossingPlane = true;
		  break;
		}
	  }
	  if(crossingPlane){break;}
	}
	for(auto& plane : twistingPlanes){
	  if(crossingPlane){break;}
	  bool firstSide = particles[c.members.front().index].position.dot(plane.normal) > plane.offset;
	  
	  for(const auto& member : c.members){
		bool thisSide = particles[member.index].position.dot(plane.normal) > plane.offset;
		if(thisSide != firstSide){
		  crossingPlane = true;
		  break;
		}
	  }
	}
	for(auto& plane : tiltingPlanes){
	  if(crossingPlane){break;}
	  bool firstSide = particles[c.members.front().index].position.dot(plane.normal) > plane.offset;
	  
	  for(const auto& member : c.members){
		bool thisSide = particles[member.index].position.dot(plane.normal) > plane.offset;
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

void World::mergeClusters(const std::vector<Particle>& newParticles,
	const std::vector<Cluster>& newClusters){

  auto particleStart = particles.size();
  auto clusterStart = clusters.size();
  
  //stick them in there and update indeces
  particles.insert(particles.end(), newParticles.begin(), newParticles.end());
  for(auto i : range(particleStart, particles.size())){
	particles[i].id = i;
	for(auto& c : particles[i].clusters){
	  c += clusterStart;
	}
  }
  
  clusters.insert(clusters.end(), newClusters.begin(), newClusters.end());
  for(auto i : range(clusterStart, clusters.size())){
	for(auto& member : clusters[i].members){
	  member.index += particleStart;
	}
  }
  
}

void World::removeInvalidClusters(){
  for(int i = 0; i < clusters.size(); i++){
    Cluster cluster = clusters.at(i);
    Cluster* toDelete = NULL;
    int index = 0;

    bool tooManyClusters = true;
    for(int j = 0; j < particles.size(); j++){
      if(particles[j].clusters.size() <= Mmax){
	tooManyClusters = false;
      }
    }

    if(cluster.members.size() < cluster.initialMembers / 2 ||
       cluster.members.size() > cluster.initialMembers * 2 ||
       cluster.Fp.norm() > 2.0 ||
       tooManyClusters){
      toDelete = &cluster;
      index = i;
    }

    if(toDelete != NULL){
      clusters.erase(clusters.begin() + index);
    }
  }
}


