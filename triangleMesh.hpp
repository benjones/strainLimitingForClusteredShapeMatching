#pragma once

#include "particle.h"
#include <array>


class TriangleMesh{
public:
  using RMMatrix3f = Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>;
  using RMMatrix3i = Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>;


  RMMatrix3f restVertices;
  RMMatrix3i triangles;

  //number of particle neighbors used to skin
  std::vector<std::array<std::pair<int, float>, 8 > > skinningWeights;
  int numParticles;

  TriangleMesh() :numParticles(0){}
  TriangleMesh(const std::string& filename) :numParticles(0){loadFromFile(filename); }

  void loadFromFile(const std::string& filename);

  void computeSkinningWeights(const std::vector<Particle>& particles);

  void deformAndOuput(const std::string& filename, const std::vector<Particle>& particles) ;
  
};
