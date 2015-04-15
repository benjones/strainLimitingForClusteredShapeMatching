
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include "json/json.h"

#include "range.hpp"
using benlib::range;


Eigen::Vector3d getSupportingPoint(const Eigen::Vector3d& normal, double offset){
  return normal*offset;
}

/*std::pair<Eigen::Vector3d, Eigen::Vector3d>
getPlaneTangents(const Eigen::Vector3d& normal){

  const Eigen::Vector3f xVector{1,0,0};
  Eigen::Vector3f span1 = (fabs(normal.dot(xVector)) < 0.9) ?
	normal.cross(xVector) :
	normal.cross(Eigen::Vector3f{0,1,0});
  Eigen::Vector3f span2 = normal.cross(span1);
  
  span1.normalize();
  span2.normalize();
  
  assert(fabs(normal.squaredNorm() - 1) < 0.01);
  assert(fabs(span1.squaredNorm() - 1) < 0.01);
  assert(fabs(span2.squaredNorm() - 1) < 0.01);
  assert(fabs(normal.dot(span1)) < .01);
  assert(fabs(normal.dot(span2)) < .01);
  assert(fabs(span1.dot(span2)) < .01);
  
  return std::make_pair(span1, span2);
}*/

int main(int argc, char** argv){


  const std::string usage = "mitsubafy <jsonFile> <formatStringForParticles> <outputDirectory> <radius/skinnedObjFormatString> [extraXMLStuff]";
  
  if(argc < 5){ std::cout << usage << std::endl; return 1;}
  
  std::string extraXML;
  if(argc == 6){
	std::ifstream ins(argv[5]);
	extraXML.assign(std::istreambuf_iterator<char>(ins),
		std::istreambuf_iterator<char>());
  }

  const std::string outputDir = argv[3];
  
  char* strEnd;
  const double radius = std::strtod(argv[4], &strEnd);
  const bool renderSpheres = (strEnd != argv[4]);
  std::string objFormatString(argv[4]);

  std::vector<std::string> planeStrings;
   
  
  //load json stuff for collision geo
  Json::Value root;
  {
	Json::Reader jReader;
	std::ifstream ins(argv[1]);
	if(!jReader.parse(ins, root)){
	  std::cout << "couldn't read input file: " << argv[1] << '\n'
				<< jReader.getFormattedErrorMessages() << std::endl;
	  return 1;
	}


	for(auto i : range(root["planes"].size())){
	  auto& plane = root["planes"][i];
	  Eigen::Vector3d normal(plane[0].asDouble(), plane[1].asDouble(), plane[2].asDouble());
	  normal.normalize();
	  double offset = plane[3].asDouble();
	  Eigen::Vector3d zAxis(0,0,1);
	  Eigen::Vector3d cp = zAxis.cross(normal);
	  
	  std::string planeString = "<shape type=\"rectangle\" >\n<transform name=\"toWorld\">\n";
	  planeString += "<scale x=\"50\" y=\"50\" z=\"50\" />\n";
	  if(cp.norm() > 1e-3){
		Eigen::Vector3d direction = cp.normalized();
		double angle = asin(std::max(-1.0, std::min(1.0, cp.norm())))*180.0/M_PI;
		planeString += "<rotate x=\"" + std::to_string(direction[0]) +
		  "\" y=\"" + std::to_string(direction[1]) +
		  "\" z=\"" + std::to_string(direction[2]) + 
		  "\" angle=\"" + std::to_string(angle) + "\" />\n";
	  }
	  auto supportingPoint = getSupportingPoint(normal, offset);
	  planeString += "<translate x=\"" + std::to_string(supportingPoint.x()) +
		"\" y=\"" + std::to_string(supportingPoint.y()) +
		"\" z=\"" + std::to_string(supportingPoint.z()) + 
		"\" />\n";
	  planeString += "</transform>\n";
	  
	  planeString += "<bsdf type=\"diffuse\"><srgb name=\"reflectance\" value=\"#444444\" /></bsdf>\n</shape>\n";
	  planeStrings.push_back(std::move(planeString));
	}






  }


  const std::string mitsubaHeader = 
	"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<scene version=\"0.5.0\">\n";
  const std::string mitsubaFooter = 
	"</scene>\n";

  const std::string shadingInfo =
	"<bsdf type=\"diffuse\"><srgb name=\"reflectance\" value=\"#aaaaaa\" /></bsdf>";

  const std::string sphereStart = "<shape type=\"sphere\">\n<point name=\"center\" ";
  const std::string sphereEnd = "</shape>\n";

  char currentFile[2048];
  for(int frame = 0; ; ++frame){
	std::vector<float> positions;
	sprintf(currentFile, argv[2], frame);
	//this will tell us how many files there are, even if we're using objs.
	std::ifstream particleIns(currentFile, std::ios_base::binary | std::ios_base::in);
	if(!particleIns.good()){
	  std::cout << "no more particle files, " << currentFile << std::endl;
	  break;
	}
	std::cout << "processing file: " << currentFile << std::endl;
	
	
	char frameString[8];
	sprintf(frameString, "%04d", frame); //why I can't use std::to_string for this is beyond me

	const std::string outputFile = outputDir + "/mitsubaFrame_" + frameString + ".xml";
	std::ofstream outs(outputFile);
	
	outs << mitsubaHeader;
	outs << extraXML;
	for(auto&& ps : planeStrings){ outs << ps << std::endl;}

	//deal with moving obstacles:
	for(auto i : range(root["movingPlanes"].size())){  

	  auto& plane = root["movingPlanes"][i];

	  Eigen::Vector3d normal(plane["normal"][0].asDouble(),
		  plane["normal"][1].asDouble(),
		  plane["normal"][2].asDouble());
	  normal.normalize();
	  double offset = plane["offset"].asDouble();
	  double velocity = plane["velocity"].asDouble();
	  Eigen::Vector3d zAxis(0,0,1);
	  Eigen::Vector3d cp = zAxis.cross(normal);
	  

	  outs << "<shape type=\"rectangle\" >\n<transform name=\"toWorld\">\n"
		   << "<scale x=\"1\" y=\"1\" z=\"1\" />\n";
	  if(cp.norm() > 1e-3){
		Eigen::Vector3d direction = cp.normalized();
		double angle = asin(std::max(-1.0, std::min(1.0, cp.norm())))*180.0/M_PI;
		outs << "<rotate x=\"" << direction[0] 
			 << "\" y=\""  << direction[1]
			 << "\" z=\""  << direction[2]  
			 << "\" angle=\""  << angle << "\" />\n";
	  }
	  auto supportingPoint = getSupportingPoint(normal, 
		  frame*root.get("dt", 1.0/60.0).asDouble()*velocity + offset);
	  outs << "<translate x=\"" << supportingPoint.x()
		   << "\" y=\""  << supportingPoint.y()
		   << "\" z=\""  << supportingPoint.z()
		   << "\" />\n"
		   << "</transform>\n"
		   << "<bsdf type=\"thindielectric\"><srgb name=\"specularTransmittance\" value=\"#ff8888\" />\n"
		   << "<srgb name=\"specularReflectance\" value=\"#333333\" />"
		   << "</bsdf>\n</shape>\n";
	

	}

	//twisting planes
	for(auto i : range(root["twistingPlanes"].size())){
	  auto& plane = root["twistingPlanes"][i];
	  std::cout << "twisting plane: " << i << std::endl;
	  double lifetime = plane.get("lifetime", 1e100).asDouble();
	  double timeNow = frame*root.get("dt", 1.0/60.0).asDouble();
	  std::cout << "lifetime: " << lifetime << " time now: " << timeNow << std::endl;
	  if(timeNow < lifetime ){
		std::cout << "active plane" << std::endl;
		Eigen::Vector3d normal(plane["normal"][0].asDouble(),
			plane["normal"][1].asDouble(),
			plane["normal"][2].asDouble());
		normal.normalize();
		double offset = plane["offset"].asDouble();
		double velocity = plane["velocity"].asDouble();
		Eigen::Vector3d zAxis(0,0,1);
		Eigen::Vector3d cp = zAxis.cross(normal);
		
		double width = plane["width"].asDouble();


		outs << "<shape type=\"rectangle\" >\n<transform name=\"toWorld\">\n"
			 << "<scale x=\"" << width << "\" y=\"" << width << "\" />\n";
		if(cp.norm() > 1e-3){
		  Eigen::Vector3d direction = cp.normalized();
		  double angle = asin(std::max(-1.0, std::min(1.0, cp.norm())))*180.0/M_PI;
		  outs << "<rotate x=\"" << direction[0] 
			   << "\" y=\""  << direction[1]
			   << "\" z=\""  << direction[2]  
			   << "\" angle=\""  << angle << "\" />\n";
		}
		//do the time dependent rotation
		double angVel = plane["angularVelocity"].asDouble();
		outs << "<rotate x=\"" << normal.x() 
			 << "\" y=\"" << normal.y()
			 << "\" z=\"" << normal.z() 
			 << "\" angle=\"" << timeNow*angVel*180/M_PI
			 << "\" />\n";
		Eigen::Vector3d supportingPoint = normal*offset;
		outs << "<translate x=\"" << supportingPoint.x()
			 << "\" y=\""  << supportingPoint.y()
			 << "\" z=\""  << supportingPoint.z()
			 << "\" />\n"
			 << "</transform>\n"
			 << "<bsdf type=\"thindielectric\"><srgb name=\"specularTransmittance\" value=\"#ff8888\" />\n"
			 << "<srgb name=\"specularReflectance\" value=\"#333333\" />"
			 << "</bsdf>\n</shape>\n";



	  }
	}


	//projectiles
	for(auto i : range(root["projectiles"].size())){
	  auto& projectile = root["projectiles"][i];
	  Eigen::Vector3d start(projectile["start"][0].asDouble(),
		  projectile["start"][1].asDouble(),
		  projectile["start"][2].asDouble());

	  Eigen::Vector3d velocity(projectile["velocity"][0].asDouble(),
		  projectile["velocity"][1].asDouble(),
		  projectile["velocity"][2].asDouble());

	  Eigen::Vector3d position = start + frame*root.get("dt", 1.0/60.0).asDouble()*velocity;

	  double sphereRadius = projectile["radius"].asDouble();
	  
	  outs << "<shape type=\"sphere\" >\n<point name=\"center\" x=\""
		   << position.x() << "\" y=\""
		   << position.y() << "\" z=\""
		   << position.z() << "\" />\n<float name=\"radius\" value=\""
		   << sphereRadius << "\" />\n"
		   << "<bsdf type=\"diffuse\"><srgb name=\"reflectance\" value=\"#3333ee\" />"
		   << "</bsdf></shape>\n";
	}


	if(renderSpheres){
	  size_t nParticles;
	  particleIns.read(reinterpret_cast<char*>(&nParticles), sizeof(nParticles));
	  positions.resize(nParticles*3);
	  particleIns.read(reinterpret_cast<char*>(positions.data()), 
		  3*nParticles*sizeof(typename decltype(positions)::value_type));
	  std::ifstream colorFile(std::string(currentFile) + ".colors", std::ios_base::binary | std::ios_base::in);
	  bool colors = colorFile.good();
	  std::vector<float> colorList;
	  if(colors){
		std::cout << "found colors" << std::endl;
		size_t numParticles;
		colorFile.read(reinterpret_cast<char*>(&numParticles), sizeof(numParticles));
		assert(numParticles == nParticles);
		colorList.resize(numParticles*3);
		colorFile.read(reinterpret_cast<char*>(colorList.data()), 
			3*numParticles*sizeof(float));
		
	  }
	  for(auto i : range(nParticles)){
		outs << sphereStart << "x=\"" 
			 << positions[3*i] << "\" y=\"" 
			 << positions[3*i + 1] << "\" z=\""
			 << positions[3*i + 2] << "\" />\n<float name=\"radius\" value=\""
			 << radius << "\" />";
		if(colors){
		  outs << "<bsdf type=\"diffuse\"><rgb name=\"reflectance\" value=\""
			   << colorList[3*i] << ", " << colorList[3*i +1] << ", " << colorList[3*i +2]
			   << "\" /></bsdf>\n";
		} else {
		  outs << shadingInfo;
		}
		outs << sphereEnd << std::endl;
	  }
	} else {
	  char objName[2048];
	  sprintf(objName, objFormatString.c_str(), frame);
	  outs << "<shape type=\"obj\" >\n"
		   << "<string name=\"filename\" value=\""
		   << objName << "\" />\n"
		   << "<bsdf type=\"diffuse\"><srgb name=\"reflectance\" value=\"#33ee33\" /></bsdf>\n"
		   << "</shape>\n";
	}
	outs << mitsubaFooter << std::endl;
  }

  return 0;
}
