
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


  const std::string usage = "mitsubafy <jsonFile> <formatStringForParticles> <outputDirectory> <radius> [extraXMLStuff]";
  
  if(argc < 5){ std::cout << usage << std::endl; return 1;}
  
  std::string extraXML;
  if(argc == 6){
	std::ifstream ins(argv[5]);
	extraXML.assign(std::istreambuf_iterator<char>(ins),
		std::istreambuf_iterator<char>());
  }

  const std::string outputDir = argv[3];
  const double radius = std::stod(argv[4]);


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
	sprintf(currentFile, argv[2], frame);
	std::ifstream particleIns(currentFile, std::ios_base::binary | std::ios_base::in);
	if(!particleIns.good()){
	  std::cout << "no more particle files, " << currentFile << std::endl;
	  break;
	}
	std::cout << "processing file: " << currentFile << std::endl;
	size_t nParticles;
	particleIns.read(reinterpret_cast<char*>(&nParticles), sizeof(nParticles));
	std::vector<float> positions(nParticles*3);
	particleIns.read(reinterpret_cast<char*>(positions.data()), 
		3*nParticles*sizeof(typename decltype(positions)::value_type));

	char frameString[8];
	sprintf(frameString, "%04d", frame); //why I can't use std::to_string for this is beyond me

	const std::string outputFile = outputDir + "/mitsubaFrame_" + frameString + ".xml";
	std::ofstream outs(outputFile);
	
	outs << mitsubaHeader;
	for(auto&& ps : planeStrings){ outs << ps << std::endl;}
	for(auto i : range(nParticles)){
	  outs << sphereStart << "x=\"" 
		   << positions[3*i] << "\" y=\"" 
		   << positions[3*i + 1] << "\" z=\""
		   << positions[3*i + 2] << "\" />\n<float name=\"radius\" value=\""
		   << radius << "\" />"
		   << shadingInfo << sphereEnd << std::endl;
	}
	outs << mitsubaFooter << std::endl;
  }

  return 0;
}
