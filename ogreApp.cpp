
#include <memory>
#include <iostream>
#include <OgreRoot.h>
#include <OgreOctreePlugin.h>
#include <OgreGLPlugin.h>
#include <OgreWindowEventUtilities.h>
#include <OgreRenderWindow.h>
#include <OgreCamera.h>
#include <OgreViewport.h>
#include <OgreSceneNode.h>
#include <OgreEntity.h>

int main(){



  auto octreePlugin = std::unique_ptr<Ogre::OctreePlugin>(new Ogre::OctreePlugin());
  auto glPlugin = std::unique_ptr<Ogre::GLPlugin>(new Ogre::GLPlugin());
  //the top two must outlive root
  auto ogreRoot = std::unique_ptr<Ogre::Root>(new Ogre::Root("","","ogreLog.log"));

  ogreRoot->installPlugin(octreePlugin.get());
  octreePlugin->initialise();


  ogreRoot->installPlugin(glPlugin.get());
  glPlugin->initialise();
  
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

  assert(rs);
  if(!(rs->getName() == "OpenGL Rendering Subsystem")){
	throw std::runtime_error("couldn't set up render sys");
  }

  rs->setConfigOption("Full Screen", "No");
  rs->setConfigOption("VSync", "No");
  rs->setConfigOption("Video Mode", "800 x 600 @ 32-bit");

  ogreRoot->setRenderSystem(rs);

  auto* window = ogreRoot->initialise(true, "Ductile Fracture for Shape Matching");

  sceneManager->setAmbientLight(Ogre::ColourValue(0.5, 0.5, 0.5));

  auto camera = sceneManager->createCamera("theCamera");
  //todo use stuff from runSimulator
  camera->setPosition(Ogre::Vector3(0,0,10));
  camera->lookAt(Ogre::Vector3(0,0,0));
  camera->setNearClipDistance(0.5);

  auto* viewport = window->addViewport(camera);
  viewport->setBackgroundColour(Ogre::ColourValue(0,0,0));

  camera->setAspectRatio(
	  Ogre::Real(viewport->getActualWidth()) /
	  Ogre::Real(viewport->getActualHeight()));

  //load the meshes/textures/etc
  Ogre::ResourceGroupManager::getSingleton().
	addResourceLocation("/Users/ben/libs/ogre/Samples/Media/models", "FileSystem", "General");
  Ogre::ResourceGroupManager::getSingleton().initialiseAllResourceGroups();

  Ogre::MovableObject* sphereEntity = sceneManager->createEntity("Sphere", "sphere.mesh");

  auto* headNode = sceneManager->getRootSceneNode()->createChildSceneNode();
  headNode->attachObject(sphereEntity);
   
  
  
  bool readyToExit = false;
  while(!readyToExit){
	// Pump window messages for nice behaviour
	Ogre::WindowEventUtilities::messagePump();
	
	// Render a frame
	ogreRoot->renderOneFrame();
	
	if(window->isClosed()){
	  readyToExit = true;
	}
  }
  
  
  return 0;
}
