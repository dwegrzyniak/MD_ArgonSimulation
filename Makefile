CC = g++
LIBS = -lncurses
GFLAGS = -Wall -pedantic -O2 -std=c++17
all:	argonSimulation
argonSimulation:	kms01.o parameters.o moleculesSystem.o molecule.o
	$(CC) -o argonSimulation kms01.o parameters.o moleculesSystem.o molecule.o $(GFLAGS) $(LIBS)
	
molecule.o: molecule.cpp molecule.hpp
	$(CC) -o molecule.o -c molecule.cpp $(GFLAGS)
	
parameters.o:	parameters.cpp molecule.hpp
	$(CC) -o parameters.o -c parameters.cpp $(GFLAGS)

moleculesSystem.o: moleculesSystem.cpp molecule.hpp
	$(CC) -o moleculesSystem.o -c moleculesSystem.cpp $(GFLAGS)

kms01.o:	kms01.cpp molecule.hpp
	$(CC) -o kms01.o -c kms01.cpp $(GFLAGS)
