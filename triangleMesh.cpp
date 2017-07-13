#include "triangleMesh.hpp"

#include <fstream>
#include "plyIO.hpp"
#include "range.hpp"
using benlib::range;

#include "accelerationGrid.h"

void TriangleMesh::loadFromFile(const std::string& filename){

  std::ifstream ins(filename, std::ios_base::binary);
  readPLY(ins, restVertices, triangles);

  
}

struct ePosGetter{
  Eigen::Vector3d operator()(const Particle& p) const{
	return p.embeddedPosition;
  }  
};

void TriangleMesh::computeSkinningWeights(const std::vector<Particle>& particles){

  numParticles = particles.size();


  
  AccelerationGrid<Particle, ePosGetter> grid;
  grid.updateGrid(particles);

  skinningWeights.resize(restVertices.rows());

  for(auto i : range(particles.size())){
	float searchRadius = 0.1;
	Eigen::Vector3d vertexPosition = restVertices.row(i).template cast<double>();
	while(true){
	  auto neighbors = grid.getNearestNeighbors(particles, vertexPosition, searchRadius);
	  if(neighbors.empty()){
		searchRadius *= 2;
	  } else {
		std::sort(neighbors.begin(), neighbors.end(),
			[&vertexPosition, &particles](int a, int b){
			  const auto& va = particles[a].embeddedPosition;
			  const auto& vb = particles[b].embeddedPosition;
			  return (va - vertexPosition).squaredNorm() <
				(vb - vertexPosition).squaredNorm();
			});
		int j = 0;
		float weightSum = 0;
		for(int j = 0; j < skinningWeights[i].size() && j < neighbors.size(); ++j){
		  float weight = std::exp(-(vertexPosition -
				  particles[neighbors[j]].embeddedPosition).squaredNorm());
		  weightSum += weight;
		  skinningWeights[i][j] = {neighbors[j], weight};
		}
		for( ; j < skinningWeights[i].size(); ++j){
		  skinningWeights[i][j] = {-1, 0};
		}
		for(auto & sw : skinningWeights[i]){
		  sw.second /= weightSum;
		}
		break;
	  }
	}
  }
}


void TriangleMesh::deformAndOuput(const std::string& filename,
	const std::vector<Particle>& particles) {

  if(particles.size() != numParticles){
	std::cout << "numParticles changed, recomputing skinning weights" << std::endl;
	computeSkinningWeights(particles);
  }
  
  RMMatrix3f deformedVertices = RMMatrix3f::Zero(restVertices.rows(), 3);
  for(auto i : range(restVertices.rows())){
	deformedVertices.row(i) = restVertices.row(i);
	for(const auto& pr : skinningWeights[i]){
	  if(pr.first != -1){
		Eigen::Vector3f diff = (particles[pr.first].position -
			particles[pr.first].embeddedPosition).template cast<float>();
		deformedVertices.row(i) += pr.second*diff;
	  }
	}
  }
  std::ofstream plyOuts(filename, std::ios_base::binary);
  writePLY(plyOuts, deformedVertices, triangles);
}
