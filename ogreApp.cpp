/*
  hg clone https://bitbucket.org/sinbad/ogre/
  CMakeLists.txt (L327): set(CMAKE_OSX_DEPLOYMENT_TARGET 10.12)
  cmake -D OGRE_STATIC=1 -D OGRE_BUILD_SAMPLES=0 ..
  make     
  make install
  cd sdk/include/OGRE
  cp OgrePlugin.h RenderSystems/GL/
  cp OgrePlugin.h Plugins/OctreeSceneManager/
  cd ../../../
  sudo mv sdk /opt/ogre
 */

#include <memory>
#include <iostream>
#ifndef BEN
#include <Ogre/OgreRoot.h>
#include <Ogre/OgrePlugin.h>
#include <Ogre/OgreWindowEventUtilities.h>
#include <Ogre/OgreRenderWindow.h>
#include <Ogre/OgreCamera.h>
#include <Ogre/OgreViewport.h>
#include <Ogre/OgreSceneNode.h>
#include <Ogre/OgreEntity.h>
#include <Ogre/OgreMaterialManager.h>
#include <Ogre/OgreMaterial.h>
#include <Ogre/OgreTechnique.h>
#include <Ogre/OgreFreeimageCodec.h>
#include <Ogre/OgreMeshManager.h>
#include <Ogre/Plugins/OctreeSceneManager/OgreOctreePlugin.h>
#include <Ogre/RenderSystems/GL/OgreGLPlugin.h>
#else
#include <Ogre/OgreRoot.h>
#include <Ogre/OgrePlugin.h>
#include <Ogre/OgreWindowEventUtilities.h>
#include <Ogre/OgreRenderWindow.h>
#include <Ogre/OgreCamera.h>
#include <Ogre/OgreViewport.h>
#include <Ogre/OgreSceneNode.h>
#include <Ogre/OgreEntity.h>
#include <Ogre/OgreMaterialManager.h>
#include <Ogre/OgreMaterial.h>
#include <Ogre/OgreTechnique.h>
#include <OgreFreeimageCodec.h>
#include <Ogre/OgreMeshManager.h>
#include <OgreOctreePlugin.h>
#include <OgreGLPlugin.h>
//#include <OgreGL3PlusPlugin.h>

#endif

#include "world.h"
#include "color_spaces.h"
#include "particle.h"
#include "range.hpp"

int main(int argc, char** argv){
  if(argc < 2){
	std::cout << "usage: ./ogreApp <input.json> [dump eye ex ey ex [lx ly lz]]" << std::endl;
	return 1;
  }
  
  
  int identity_id=0;

  bool dumpFrames = argc > 2;

  Eigen::Vector3d eye{0,0,3};
  if(argc >= 7){
	eye.x() = atof(argv[4]);
	eye.y() = atof(argv[5]);
	eye.z() = atof(argv[6]);
  }

  Eigen::Vector3d lPos{20,25,20};
  if(argc >= 10){
	lPos.x() = atof(argv[7]);
	lPos.y() = atof(argv[8]);
	lPos.z() = atof(argv[9]);
  }
  

  auto octreePlugin = std::unique_ptr<Ogre::OctreePlugin>(new Ogre::OctreePlugin());
  auto glPlugin = std::unique_ptr<Ogre::GLPlugin>(new Ogre::GLPlugin());
  //  auto gl3Plugin = std::unique_ptr<Ogre::GL3PlusPlugin>(new Ogre::GL3PlusPlugin());
  //the top two must outlive root
  auto ogreRoot = std::unique_ptr<Ogre::Root>(new Ogre::Root("","","ogreLog.log"));

  Ogre::FreeImageCodec::startup();

  
  ogreRoot->installPlugin(octreePlugin.get());
  octreePlugin->initialise();

  ogreRoot->installPlugin(glPlugin.get());
  glPlugin->initialise();
  //  ogreRoot->installPlugin(gl3Plugin.get());
  //  gl3Plugin->initialise();
  
  std::cout << "installed plugins: " << std::endl;
  for(const auto& plugin : ogreRoot->getInstalledPlugins()){
	std::cout << '\t' << plugin->getName() << std::endl;
  }

  auto* sceneManager = ogreRoot->createSceneManager("OctreeSceneManager");

  auto renderers = ogreRoot->getAvailableRenderers();
  std::cout << "renderers: " << std::endl;
  for(auto & r : renderers){
	std::cout << '\t' << r->getName() << std::endl;
  }
  
Ogre::RenderSystem* rs =
										 ogreRoot->getRenderSystemByName("OpenGL Rendering Subsystem");
	//	ogreRoot->getRenderSystemByName("OpenGL 3+ Rendering Subsystem");

  assert(rs);
  if(!(rs->getName() == "OpenGL Rendering Subsystem")){
	throw std::runtime_error("couldn't set up render sys");
  }
std::cout << "setting fullscreen" << std::endl;
  rs->setConfigOption("Full Screen", "No");
  //rs->setConfigOption("VSync", "No");
std::cout << "setting video mode" << std::endl;
  rs->setConfigOption("Video Mode", "960 x 540 @ 32-bit");
std::cout << "setting render system: " << std::endl;

ogreRoot->setRenderSystem(rs);
std::cout << "initializing window " << std::endl;

auto* window = ogreRoot->initialise(true, "Ductile Fracture for Shape Matching");

  sceneManager->setAmbientLight(Ogre::ColourValue(0.1, 0.1, 0.1));

  auto camera = sceneManager->createCamera("theCamera");
  //todo use stuff from runSimulator
  camera->setPosition(Ogre::Vector3(eye.x(),eye.y(),eye.z()));
  camera->lookAt(Ogre::Vector3(0,0,0));//Ogre::Vector3(0,-2,0));
  camera->setNearClipDistance(0.5);

  auto* viewport = window->addViewport(camera);
  viewport->setBackgroundColour(Ogre::ColourValue(0.3,0.3,0.3));

  camera->setAspectRatio(
	  Ogre::Real(viewport->getActualWidth()) /
	  Ogre::Real(viewport->getActualHeight()));

  //load the meshes/textures/etc
  Ogre::ResourceGroupManager::getSingleton().
	addResourceLocation("/opt/ogre/Media/models", "FileSystem", "General");
  Ogre::ResourceGroupManager::getSingleton().initialiseAllResourceGroups();


  auto* keyLight = sceneManager->createLight("key");

  
  
  keyLight->setPosition(lPos.x(), lPos.y(), lPos.z());


  //load stuff
  World world;
  world.loadFromJson(argv[1]);
  world.initializeNeighbors();

  const double meshSize = 100;
  const double sphereSize = 0.01; //.03 for broken heart, armadillo
  const double scaleFactor = sphereSize/meshSize;

  //SM should outlive particles...
  auto cleanupParticle =
	[&sceneManager](Particle& p){
	/*if(p.sceneNode != nullptr){
	  std::cout << "make invisible" << std::endl;
	  p.sceneNode->setVisible(false);
	}
	return;*/
	if(p.entity == nullptr){return;}
	p.sceneNode->detachObject(p.entity);
	sceneManager->destroySceneNode(p.sceneNode);
	sceneManager->destroyEntity(p.entity);
	//std::cout << "cleaning up a particle " << std::endl;
	
	
  };


  auto projectileMaterial = Ogre::MaterialManager::getSingleton().create(
	  "projMat",
	  Ogre::ResourceGroupManager::DEFAULT_RESOURCE_GROUP_NAME);

  projectileMaterial->getTechnique(0)->setDiffuse(0,0,1,1);
  
  
  std::vector<std::pair<Ogre::Entity*, Ogre::SceneNode*> > ogreProjectiles;
  for(auto & p : world.projectiles){
	ogreProjectiles.emplace_back();
	auto& op = ogreProjectiles.back();
	op.second = sceneManager->getRootSceneNode()->createChildSceneNode();
	op.second->setScale(p.radius/meshSize,p.radius/meshSize,p.radius/meshSize);
	op.second->setPosition(p.start.x(), p.start.y(), p.start.z());
	op.first = sceneManager->createEntity("sphere.mesh");
	op.second->attachObject(op.first);
	op.first->setMaterial(projectileMaterial);
  }

  auto planeMaterial = Ogre::MaterialManager::getSingleton().create(
	  "planeMat",
	  Ogre::ResourceGroupManager::DEFAULT_RESOURCE_GROUP_NAME);
  planeMaterial->getTechnique(0)->setDiffuse(1,0,0,1);
  

  std::vector<std::pair<Ogre::Entity*, Ogre::SceneNode*> > ogreMovingPlanes;
  int planeCount = 0;
  for(auto & mp : world.movingPlanes){
	
	std::string name = "movingPlane" + std::to_string(planeCount++);
	std::cout << "making moving plane: " << name << std::endl;
	ogreMovingPlanes.emplace_back();
	auto& omp = ogreMovingPlanes.back();
	omp.second = sceneManager->getRootSceneNode()->createChildSceneNode();
	Ogre::Plane plane(-Ogre::Vector3(mp.normal.x(), mp.normal.y(), mp.normal.z()), -mp.offset);
	if(fabs(1.0 - fabs(mp.normal.dot(Eigen::Vector3d(0,1,0)))) < .001){
	  std::cout << "pointing up" << std::endl;
	  //it points up, ogre is dumb here...
	  Ogre::MeshManager::getSingleton().createPlane(
		  name,
		  Ogre::ResourceGroupManager::DEFAULT_RESOURCE_GROUP_NAME,
		  plane,
		  2, 2,
		  20,20,
		  true, 1, 5, 5,
		  Ogre::Vector3::UNIT_Z);


	} else {
	  std::cout << "pointing somewhere else" << std::endl;
	  Ogre::MeshManager::getSingleton().createPlane(
		  name,
		  Ogre::ResourceGroupManager::DEFAULT_RESOURCE_GROUP_NAME,
		  plane,
		  2, 2,
		  20,20);
	}
	omp.first = sceneManager->createEntity(name);
	omp.second->attachObject(omp.first);
	omp.first->setMaterial(planeMaterial);
	

	
  }

  planeCount = 0;
  for(auto& plane : world.planes){

  }
  
  
  
  
  for(auto& p : world.particles){
	p.cleanup = cleanupParticle;
	p.sceneNode = sceneManager->getRootSceneNode()->createChildSceneNode();
	p.sceneNode->setScale(scaleFactor, scaleFactor, scaleFactor);
	p.sceneNode->setPosition(p.position.x(), p.position.y(), p.position.z());
	p.entity = sceneManager->createEntity("sphere.mesh");
	p.sceneNode->attachObject(p.entity);
	auto material = Ogre::MaterialManager::getSingleton().create(
		(std::string("aMat ") + std::to_string(identity_id++)).c_str(),
		Ogre::ResourceGroupManager::DEFAULT_RESOURCE_GROUP_NAME);

	//find closest color:
	/*	auto closestCluster = std::min_element(
		p.clusters.begin(), p.clusters.end(),
		[&p,&world](int a, int b){
		  return (p.position - world.clusters[a].worldCom).squaredNorm() <
		  (p.position - world.clusters[b].worldCom).squaredNorm();
		  });
	
		RGBColor rgb = HSLColor(2.0*acos(-1)*(*closestCluster%12)/12.0, 0.7, 0.7).to_rgb();*/
	material->getTechnique(0)->setDiffuse(p.color.r, p.color.g, p.color.b, 1);
	p.entity->setMaterial(material);
  }
  
  bool readyToExit = false;

  std::string fileBase = "ogreFrames/frame.%04d.png";
  char fname[1024];

  int frame = 0;
  while(!readyToExit){

	world.timestep();
	for(auto& p : world.particles){
	  if(p.sceneNode == nullptr){
		p.cleanup = cleanupParticle;
		p.sceneNode = sceneManager->getRootSceneNode()->createChildSceneNode();
		p.sceneNode->setScale(scaleFactor, scaleFactor, scaleFactor);
		p.entity = sceneManager->createEntity("sphere.mesh");
		p.sceneNode->attachObject(p.entity);
		auto material = Ogre::MaterialManager::getSingleton().create(
			(std::string("aMat ") + std::to_string(identity_id++)).c_str(),
			Ogre::ResourceGroupManager::DEFAULT_RESOURCE_GROUP_NAME);
		
		//find closest color:
		/*		auto closestCluster = std::min_element(
			p.clusters.begin(), p.clusters.end(),
			[&p,&world](int a, int b){
			  return (p.position - world.clusters[a].worldCom).squaredNorm() <
			  (p.position - world.clusters[b].worldCom).squaredNorm();
			  });
		
			  RGBColor rgb = HSLColor(2.0*acos(-1)*(*closestCluster%12)/12.0, 0.7, 0.7).to_rgb();*/
		material->getTechnique(0)->setDiffuse(p.color.r, p.color.g, p.color.b, 1);
		p.entity->setMaterial(material);
	  }
	  
	  p.sceneNode->setPosition(p.position.x(), p.position.y(), p.position.z());
	  //p.sceneNode->setVisible(frame % 2 == 0);
	}

	for(auto i : benlib::range(world.projectiles.size())){
	  Eigen::Vector3d posNow = world.projectiles[i].start + world.elapsedTime*
		world.projectiles[i].velocity;
	  ogreProjectiles[i].second->setPosition(posNow.x(), posNow.y(), posNow.z());
	  
	}
	
	for(auto i : benlib::range(world.movingPlanes.size())){
	  Eigen::Vector3d positionNow = world.movingPlanes[i].velocity*world.elapsedTime*
		world.movingPlanes[i].normal;
	  ogreMovingPlanes[i].second->setPosition(positionNow.x(), positionNow.y(), positionNow.z());
	}
	

	
	// Pump window messages for nice behaviour
	Ogre::WindowEventUtilities::messagePump();
	
	// Render a frame
	ogreRoot->renderOneFrame();
	if(dumpFrames){
	  sprintf(fname, fileBase.c_str(), frame);

	  window->writeContentsToFile(fname);
	}
	std::cout << "finished frame: " << frame << std::endl;
	frame++;
	
	if(window->isClosed()){
	  readyToExit = true;
	}
  }
  
  
  return 0;
}
