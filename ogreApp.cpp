
#include <memory>
#include <OgreRoot.h>


int main(){

  auto ogreRoot = std::unique_ptr<Ogre::Root>(new Ogre::Root("","","ogreLog.log"));


  return 0;
}
