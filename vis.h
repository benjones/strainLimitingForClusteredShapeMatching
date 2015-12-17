#pragma once

#include <vector>
#include <Eigen/Dense>
#include <SDL.h>

class World;
#include "planes.hpp"

struct Camera{
  Eigen::Vector3d position = Eigen::Vector3d{0,7, -5}, 
	lookAt = Eigen::Vector3d::Zero(), 
	up = Eigen::Vector3d{0,1,0};

  void zoom(int amount);
  void pan(Eigen::Vector2i oldPosition, Eigen::Vector2i newPosition);
  void move(bool forward);
};

struct VisSettings{
  bool drawClusters = true;
  bool drawFracturePlanes = true;
  bool drawColoredParticles = false;
  bool joshDebugFlag = true;
  bool colorByToughness = false;
  int which_cluster = -1;
};

void drawWorldPretty(const World& world,
					   const Camera& camera, 
					   const VisSettings& settings, 
					   SDL_Window* window);

void drawPlanesPretty(const std::vector<Plane>& planes,
						const std::vector<MovingPlane>& movingPlanes,
						const std::vector<TwistingPlane>& twistingPlanes,
						const std::vector<TiltingPlane>& tiltingPlanes,
						double elapsedTime);

void drawTPlane(const Eigen::Vector3d& normal, double offset, double roffset, double width);
void drawTiltPlane(const Eigen::Vector3d& normal, 
	const Eigen::Vector3d& tilt, 
	double offset, double roffset, double width);


