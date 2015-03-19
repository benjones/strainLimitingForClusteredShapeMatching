CC          = clang++
cc          = clang++


FLAGS=-Wall -Wno-c++11-extensions -std=c++11 -Wno-deprecated-declarations 

#-----------------------------------------
#Optimization ----------------------------
OPT = -O2 -g

#-----------------------------------------
#-----------------------------------------

TARGETS = runSimulator

OBJECTS =  main.o particle.o world.o jsoncpp.o

#-----------------------------------------
# Update to point to your eigen headers, sdl headers
INCS = -I/Library/Frameworks/SDL2.framework/headers
EIGEN_INCLUDE=-I/usr/local/include/eigen3

# Update to point to your SDL and opengl libs
SDL_LIB = -framework SDL2 -framework OpenGL -F/Library/Frameworks



CCOPTS = $(OPT) $(FLAGS) $(INCS) $(EIGEN_INCLUDE) 
LDOPTS = $(OPT) 

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS) 


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
