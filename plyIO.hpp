#pragma once

#include "Eigen/Dense"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>



inline void writePlyHeader(std::ofstream& outs, size_t numVertices, size_t numFaces){

  outs << "ply\n"
	"format binary_little_endian 1.0\n"
	"element vertex " << numVertices
	   << "\n"
	"property float32 x\n"
	"property float32 y\n"
	"property float32 z\n"
	"element face " << numFaces
	   << "\n" 
	"property list uchar int32 vertex_indices\n"
	"end_header\n";
  

}

inline void writePlyHeaderWithTexture(std::ofstream& outs, size_t numVertices, size_t numFaces){

  outs << "ply\n"
	"format binary_little_endian 1.0\n"
	"element vertex " << numVertices
	   << "\n"
	"property float32 x\n"
	"property float32 y\n"
	"property float32 z\n"
	"property float32 u\n"
	"property float32 v\n"
	"element face " << numFaces
	   << "\n" 
	"property list uchar int32 vertex_indices\n"
	"end_header\n";
  

}


inline void writePlyHeaderWithColor(std::ofstream& outs, size_t numVertices, size_t numFaces){

  outs << "ply\n"
	"format binary_little_endian 1.0\n"
	"element vertex " << numVertices
	   << "\n"
	"property float32 x\n"
	"property float32 y\n"
	"property float32 z\n"
	"property float32 red\n"
	"property float32 green\n"
	"property float32 blue\n"
	"element face " << numFaces
	   << "\n" 
	"property list uchar int32 vertex_indices\n"
	"end_header\n";
  

}

//should work for eigen type stuff
//has to be ROW MAJOR FLOATX3 Matrix
template<typename VType, typename TType>
void writePLY(std::ofstream& outs, const VType& vertices, const TType& triangles){

  assert(vertices.cols() == 3);
  assert(triangles.cols() == 3);
  writePlyHeader(outs, vertices.rows(), triangles.rows());
  
  outs.write(reinterpret_cast<const char*>(vertices.data()),
			 3*vertices.rows()*sizeof(float));
  const unsigned char three{3u};
  for(size_t i = 0; i < triangles.rows(); ++i){
	outs.write(reinterpret_cast<const char*>(&three), sizeof(three));
	outs.write(reinterpret_cast<const char*>(triangles.data() + 3*i),
			   3*sizeof(int));
  }
}

//should work for eigen type stuff
//has to be ROW MAJOR FLOATX3 Matrix
template<typename VType, typename TType, typename TexUType, typename TexVType>
void writePLYWithTexture(std::ofstream& outs, const VType& vertices, 
	const TexUType& texU, const TexVType& texV,
	const TType& triangles){

  static_assert(std::is_same<float, typename VType::Scalar>::value, 
	  "vertices must be float type");
  /*  static_assert(std::is_same<float, typename TexUType::Scalar>::value, 
	  "textures must be float type");
  static_assert(std::is_same<float, typename TexVType::Scalar>::value, 
	  "textures must be float type");
  */
  assert(vertices.cols() == 3);
  assert(triangles.cols() == 3);
  writePlyHeaderWithTexture(outs, vertices.rows(), triangles.rows());
  
  Eigen::Matrix<float, Eigen::Dynamic, 5, Eigen::RowMajor> vertexMatrix;
  vertexMatrix.resize(vertices.rows(), 5);
  vertexMatrix.block(0, 0, vertices.rows(), 3) = vertices;
  vertexMatrix.col(3) = texU;
  vertexMatrix.col(4) = texV;

  outs.write(reinterpret_cast<const char*>(vertexMatrix.data()),
			 5*vertices.rows()*sizeof(float));

  const unsigned char three{3u};
  for(size_t i = 0; i < triangles.rows(); ++i){
	outs.write(reinterpret_cast<const char*>(&three), sizeof(three));
	outs.write(reinterpret_cast<const char*>(triangles.data() + 3*i),
			   3*sizeof(int));
  }
}

//should work for eigen type stuff
//has to be ROW MAJOR FLOATX3 Matrix
template<typename VType, typename TType, typename CType>
void writePLYWithColor(std::ofstream& outs, const VType& vertices, 
	const CType& colors,
	const TType& triangles){

  static_assert(std::is_same<float, typename VType::Scalar>::value, 
	  "vertices must be float type");
  /*  static_assert(std::is_same<float, typename TexUType::Scalar>::value, 
	  "textures must be float type");
  static_assert(std::is_same<float, typename TexVType::Scalar>::value, 
	  "textures must be float type");
  */
  assert(vertices.cols() == 3);
  assert(triangles.cols() == 3);
  writePlyHeaderWithColor(outs, vertices.rows(), triangles.rows());
  
  Eigen::Matrix<float, Eigen::Dynamic, 6, Eigen::RowMajor> vertexMatrix;
  vertexMatrix.resize(vertices.rows(), 6);
  vertexMatrix.block(0, 0, vertices.rows(), 3) = vertices;
  vertexMatrix.block(0, 3, vertices.rows(), 3) = colors;

  outs.write(reinterpret_cast<const char*>(vertexMatrix.data()),
			 6*vertices.rows()*sizeof(float));

  const unsigned char three{3u};
  for(size_t i = 0; i < triangles.rows(); ++i){
	outs.write(reinterpret_cast<const char*>(&three), sizeof(three));
	outs.write(reinterpret_cast<const char*>(triangles.data() + 3*i),
			   3*sizeof(int));
  }
}





inline void writePLY(std::ofstream& outs, 
			  const std::vector<float>& vertices, 
			  const std::vector<int>& triangles){

  assert((vertices.size() % 3) == 0);
  assert((triangles.size() % 3) == 0);

  writePlyHeader(outs, vertices.size()/3, triangles.size()/3);
  
  outs.write(reinterpret_cast<const char*>(vertices.data()),
			 vertices.size()*sizeof(float));

  const unsigned char three{3u};
  for(size_t i = 0; i < triangles.size()/3; ++i){
	outs.write(reinterpret_cast<const char*>(&three), sizeof(three));
	outs.write(reinterpret_cast<const char*>(triangles.data()+3*i), 
			   3*sizeof(int));
  }

}

//returns hasTexCoords, hasColors, numVertices, numFaces

std::tuple<bool, bool, int, int> 
inline readPLYHeader(std::ifstream& ins){
  std::string magic;
  ins >> magic;
  assert(magic == "ply");
  std::string firstWord;
  //read until we get to the end of the header
  int numVertices{-1};
  int numFaces{-1};
  bool hasTextures = false;
  bool hasColors = false;
  for(std::string line; std::getline(ins, line);){
	//std::cout << "read line: " << line << std::endl;
	std::stringstream lineStream(line);
	lineStream >> firstWord;
	if(firstWord == std::string("element")){
	  std::string secondWord;
	  lineStream >> secondWord;
	  //std::cout << "second word: " << secondWord << std::endl;
	  if(secondWord == "vertex"){
		lineStream >> numVertices;
	  } else if(secondWord == "face"){
		lineStream >> numFaces;
	  }
	} else if(firstWord == std::string("property")){
	  std::string secondWord, thirdWord;
	  lineStream >> secondWord >> thirdWord;
	  if(thirdWord == "u" || thirdWord == "v"){
		hasTextures = true;
	  } else if(thirdWord == "red" || thirdWord == "green" || thirdWord == "blue"){
		hasColors = true;
	  }
	  
	}else if(firstWord == std::string("end_header")){
	  break;
	}
	//ignore anythign else
	
  }
  std::cout << "has texture: " << hasTextures << std::endl;
  return std::make_tuple(hasTextures, hasColors, numVertices, numFaces);
}




/*this will work for files that I wrote, but probably not much else.
Assumes only 3D vertices with floats, and triangles, with int indices
and nothing else
MATRICES NEED TO BE ROW MAJOR!!!
*/
template<typename VType, typename TType>
void readPLY(std::ifstream& ins, VType& vertices, TType& triangles){
  static_assert(std::is_same<float, typename VType::Scalar>::value, 
	  "Vertices must be type float");
  static_assert(std::is_same<int,   typename TType::Scalar>::value,
	  "Triangles must be type int");


  int numVertices, numFaces;
  bool hasTexture, hasColors;
  std::tie(hasTexture, hasColors, numVertices, numFaces) = readPLYHeader(ins);
  if(numVertices <= 0){ 
	vertices.resize(0,3);
	triangles.resize(0,3);
	return;
  }
  assert(numFaces > 0);
  vertices.resize(numVertices, 3);
  triangles.resize(numFaces, 3);

  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> vertexMatrix;
  
  int nCols = 3 + (hasColors ? 3 : 0) + (hasTexture ? 2 : 0);
  vertexMatrix.resize(numVertices, nCols);
  ins.read(reinterpret_cast<char*>(vertexMatrix.data()),
	  numVertices*nCols*sizeof(float));
  vertices = vertexMatrix.block(0,0,vertices.rows(), 3);

  unsigned char three;
  for(int i = 0; i < numFaces; ++i){
	ins.read(reinterpret_cast<char*>(&three), sizeof(three));
	assert(three == 3);
	ins.read(reinterpret_cast<char*>(triangles.data() + 3*i),
			 3*sizeof(int));
  }
  

}


/*this will work for files that I wrote, but probably not much else.
Assumes only 3D vertices with floats, and triangles, with int indices
and nothing else
MATRICES NEED TO BE ROW MAJOR!!!
*/
template<typename VType, typename TType, typename TexUType, typename TexVType>
void readPLYWithTexture(std::ifstream& ins, VType& vertices, 
	TexUType& texU, TexVType& texV,
	TType& triangles){
  static_assert(std::is_same<float, typename VType::Scalar>::value, 
	  "Vertices must be type float");
  static_assert(std::is_same<int,   typename TType::Scalar>::value,
	  "Triangles must be type int");


  int numVertices, numFaces;
  bool hasTexture, hasColors;
  std::tie(hasTexture, hasColors, numVertices, numFaces) = readPLYHeader(ins);
  assert(hasTexture);
  assert(!hasColors);
  if(numVertices <= 0){ 
	vertices.resize(0,3);
	triangles.resize(0,3);
	return;
  }
  assert(numFaces > 0);
  vertices.resize(numVertices, 3);
  triangles.resize(numFaces, 3);
  texU.resize(numVertices);
  texV.resize(numVertices);

  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> vertexMatrix;
  
  vertexMatrix.resize(numVertices, 5);
  ins.read(reinterpret_cast<char*>(vertexMatrix.data()),
	  numVertices*5*sizeof(float));
  vertices = vertexMatrix.block(0,0,vertices.rows(), 3);
  texU = vertexMatrix.col(3);
  texV = vertexMatrix.col(4);


  unsigned char three;
  for(int i = 0; i < numFaces; ++i){
	ins.read(reinterpret_cast<char*>(&three), sizeof(three));
	assert(three == 3);
	ins.read(reinterpret_cast<char*>(triangles.data() + 3*i),
			 3*sizeof(int));
  }
  

}

/*this will work for files that I wrote, but probably not much else.
Assumes only 3D vertices with floats, and triangles, with int indices
and nothing else
MATRICES NEED TO BE ROW MAJOR!!!
*/
template<typename VType, typename TType, typename CType>
void readPLYWithColor(std::ifstream& ins, VType& vertices, 
	CType& colors,
	TType& triangles){
  static_assert(std::is_same<float, typename VType::Scalar>::value, 
	  "Vertices must be type float");
  static_assert(std::is_same<int,   typename TType::Scalar>::value,
	  "Triangles must be type int");


  int numVertices, numFaces;
  bool hasTexture, hasColors;
  std::tie(hasTexture, hasColors, numVertices, numFaces) = readPLYHeader(ins);
  assert(hasColors);
  assert(!hasTexture);
  if(numVertices <= 0){ 
	vertices.resize(0,3);
	triangles.resize(0,3);
	return;
  }
  assert(numFaces > 0);
  vertices.resize(numVertices, 3);
  triangles.resize(numFaces, 3);
  colors.resize(numVertices, 3);

  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> vertexMatrix;
  
  vertexMatrix.resize(numVertices, 6);
  ins.read(reinterpret_cast<char*>(vertexMatrix.data()),
	  numVertices*6*sizeof(float));
  vertices = vertexMatrix.block(0,0,vertices.rows(), 3);
  colors = vertexMatrix.block(0,3,vertices.rows(), 3);


  unsigned char three;
  for(int i = 0; i < numFaces; ++i){
	ins.read(reinterpret_cast<char*>(&three), sizeof(three));
	assert(three == 3);
	ins.read(reinterpret_cast<char*>(triangles.data() + 3*i),
			 3*sizeof(int));
  }
  

}
