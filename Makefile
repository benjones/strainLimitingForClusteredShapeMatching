CC = clang++
cc = clang++

OS := $(shell uname)
USER := $(shell whoami)

ifeq ($(OS), Darwin)
# Run MacOS commands 
FLAGS=-Wall -Wno-c++11-extensions -std=c++11 -Wno-deprecated-declarations 
# Update to point to your eigen headers, sdl headers
INCS = -I/usr/local/opt/sdl2/include/SDL2/ -I/Library/Frameworks/SDL2.framework/Headers
EIGEN_INCLUDE=-I/usr/local/opt/eigen/include/eigen3/ -I/usr/local/include/eigen3
# Update to point to your SDL and opengl libs
SDL_LIB = -L/usr/local/opt/sdl2/lib  -framework OpenGL -F/Library/Frameworks -framework SDL2
else
# check for Linux and run other commands
FLAGS = -Wall -Wno-c++11-extensions -std=c++11 -Wno-deprecated-declarations -lGL -lSDL2 -lGLU
# Update to point to your eigen headers, sdl headers
INCS = -I /usr/include/SDL2 -I /usr/include/ 
EIGEN_INCLUDE = -I /usr/include/eigen3
# Update to point to your SDL and opengl libs
SDL_LIB = -L/usr/local/opt/sdl2/lib  #-framework OpenGL #-F/Library/Frameworks -framework SDL2
endif


#-----------------------------------------
#Optimization ----------------------------
OPT = -O3 -g

#-----------------------------------------
#-----------------------------------------

TARGETS = runSimulator mitsubafy mitsubafyClusters

OBJECTS =  particle.o world.o jsoncpp.o movingPlane.o twistingPlane.o tiltingPlane.o constraintPlane.o color_spaces.o projectile.o cylinderObstacle.o dynamics.o clustering.o

#-----------------------------------------

OGRE_INCS = -I/usr/local/include -I/opt/ogre/include
OGRE_LIBS = -framework CoreFoundation -framework Cocoa -framework OpenGL -framework AGL -L/usr/local/lib -ltbb -lfreeimage -lboost_system -L/opt/ogre/lib -lOgreMainStatic -lRenderSystem_GLStatic -lPlugin_OctreeSceneManagerStatic -lOgreGLSupportStatic -lzzip -lboost_thread-mt -lz

ifeq ($(USER), ben)
EIGEN_INCLUDE=-I/Users/ben/libs/eigen/
OGRE_INCS = -I/Users/ben/libs/ogre/OgreMain/include -I/Users/ben/libs/ogre/build/include -I/usr/local/include -I/Users/ben/libs/ogre/PlugIns/OctreeSceneManager/include -I/Users/ben/libs/ogre/RenderSystems/GL/include #-I/Users/ben/libs/ogre/RenderSystems/GL3Plus/include
OGRE_LIBS = -framework CoreFoundation -framework Cocoa -framework OpenGL -framework AGL -L/usr/local/lib -ltbb -lfreeimage -lzzip -L/Users/ben/libs/ogre/build/lib/macosx -lOgreMainStatic -lOgreGLSupportStatic -lPlugin_OctreeSceneManagerStatic -lRenderSystem_GLStatic #-lRenderSystem_GL3PlusStatic #-F/Library/Frameworks -framework Ogre
FLAGS += -DBEN
endif



CCOPTS = $(OPT) $(FLAGS) $(INCS) $(EIGEN_INCLUDE) 
LDOPTS = $(OPT) 

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS) 

mitsubafy: mitsubafy.cpp jsoncpp.o
	$(CC) $(OPT) $(FLAGS) $(EIGEN_INCLUDE) -o mitsubafy mitsubafy.cpp jsoncpp.o

mitsubafyClusters: mitsubafyClusters.cpp jsoncpp.o
	$(CC) $(OPT) $(FLAGS) $(EIGEN_INCLUDE) -o mitsubafyClusters mitsubafyClusters.cpp jsoncpp.o color_spaces.o


libshapematch.a: $(OBJECTS)
	ar rcs libshapematch.a $(OBJECTS)

clean:
	/bin/rm -fv *.o $(TARGETS) 

#-----------------------------------------
#-----------------------------------------

runSimulator: $(OBJECTS) libshapematch.a main.o vis.o
	$(CC) $(LDOPTS) -o runSimulator -L. -lshapematch $(SDL_LIB) $(EIGEN_INCLUDE) $(INCS) $(FLAGS) main.o vis.o $(OBJECTS)

ogreApp.o: ogreApp.cpp libshapematch.a
	$(CC) -c -o ogreApp.o -std=c++ $(OGRE_INCS) $(EIGEN_INCLUDE) $(FLAGS) $(INCS) -O3 -g ogreApp.cpp

ogreApp: libshapematch.a ogreApp.o
	$(CC) -o ogreApp $(OGRE_LIBS) -L. -lshapematch ogreApp.o

convergenceTest: libshapematch.a convergenceTestMain.o vis.o
	$(CC) $(LDOPTS) -o convergenceTest -L. -lshapematch $(SDL_LIB) $(EIGEN_INCLUDE) $(INCS) $(FLAGS) convergenceTestMain.o vis.o $(OBJECTS)

#-----------------------------------------

.cpp.o:
	$(CC) $(CCOPTS) -c $< -o $@

#-----------------------------------------
#-----------------------------------------
