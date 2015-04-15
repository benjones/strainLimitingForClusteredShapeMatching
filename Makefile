CC          = clang++
cc          = clang++


FLAGS=-Wall -Wno-c++11-extensions -std=c++11 -Wno-deprecated-declarations 

#-----------------------------------------
#Optimization ----------------------------
OPT = -O2 -g

#-----------------------------------------
#-----------------------------------------

TARGETS = runSimulator mitsubafy

OBJECTS =  main.o particle.o world.o jsoncpp.o movingPlane.o twistingPlane.o tiltingPlane.o color_spaces.o projectile.o cylinderObstacle.o

#-----------------------------------------
# Update to point to your eigen headers, sdl headers
INCS = -I/usr/local/opt/sdl2/include/SDL2/ -I/Library/Frameworks/SDL2.framework/headers
EIGEN_INCLUDE=-I/usr/local/opt/eigen/include/eigen3/ -I/usr/local/include/eigen3/

# Update to point to your SDL and opengl libs
SDL_LIB = -framework OpenGL -framework SDL2 -F/Library/Frameworks #-L/usr/local/opt/sdl2/lib -lSDL2 



CCOPTS = $(OPT) $(FLAGS) $(INCS) $(EIGEN_INCLUDE) 
LDOPTS = $(OPT) 

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS) 

mitsubafy: mitsubafy.cpp jsoncpp.o
	$(CC) $(OPT) $(FLAGS) $(EIGEN_INCLUDE) -o mitsubafy mitsubafy.cpp jsoncpp.o

clean:
	/bin/rm -fv *.o $(TARGETS) 

#-----------------------------------------
#-----------------------------------------

runSimulator: $(OBJECTS)
	$(CC) $(LDOPTS) -o runSimulator $(OBJECTS) $(SDL_LIB)


#-----------------------------------------

.cpp.o:
	$(CC) $(CCOPTS) -c $< -o $@

#-----------------------------------------
#-----------------------------------------
