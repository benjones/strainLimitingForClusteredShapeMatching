CC          = clang++
cc          = clang++


FLAGS=-Wall -Wno-c++11-extensions -std=c++11 -Wno-deprecated-declarations 

#-----------------------------------------
#Optimization ----------------------------
OPT = -O2 -g

#-----------------------------------------
#-----------------------------------------

TARGETS = runSimulator mitsubafy mitsubafyClusters

OBJECTS =  particle.o world.o jsoncpp.o movingPlane.o twistingPlane.o tiltingPlane.o color_spaces.o projectile.o cylinderObstacle.o dynamics.o clustering.o

#-----------------------------------------
# Update to point to your eigen headers, sdl headers
INCS = -I/usr/local/opt/sdl2/include/SDL2/ -I/Library/Frameworks/SDL2.framework/Headers
EIGEN_INCLUDE=-I/usr/local/opt/eigen/include/eigen3/ -I/usr/local/include/eigen3

# Update to point to your SDL and opengl libs
SDL_LIB = -L/usr/local/opt/sdl2/lib  -framework OpenGL -F/Library/Frameworks -framework SDL2


OGRE_INCS = -I/Users/ben/libs/ogre/OgreMain/include -I/Users/ben/libs/ogre/build/include -I/usr/local/include
OGRE_LIBS = -framework CoreFoundation -framework Cocoa -L/usr/local/lib -ltbb -L/Users/ben/libs/ogre/build/lib/macosx -lOgreMainStatic #-F/Library/Frameworks -framework Ogre



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
	$(CC) $(LDOPTS) -o runSimulator -L. -lshapematch $(SDL_LIB) $(EIGEN_INCLUDE) $(INCS) $(FLAGS) main.o vis.o

ogreApp.o: ogreApp.cpp libshapematch.a
	$(CC) -c -o ogreApp.o -std=c++ $(OGRE_INCS) $(FLAGS) $(OPT) ogreApp.cpp

ogreApp: libshapematch.a ogreApp.o
	$(CC) -o ogreApp $(OGRE_LIBS) ogreApp.o


#-----------------------------------------

.cpp.o:
	$(CC) $(CCOPTS) -c $< -o $@

#-----------------------------------------
#-----------------------------------------
