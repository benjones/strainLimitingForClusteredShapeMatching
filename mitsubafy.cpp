
#include <fstream>
#include <iostream>

#include "json/json.h"

#include "range.hpp"
using benlib::range;



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
