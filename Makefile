all: compile link

compile:
	g++ -c main.cpp
	g++ -c particle.cpp
	g++ -c particle_sim.cpp

# add -mwindows at end of link to hide console
link:
	g++ main.o particle.o particle_sim.o -o main -lsfml-graphics -lsfml-window -lsfml-system	
