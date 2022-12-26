//parameters: n, a, f, R, L, m
#include "parameters.hpp"
#include "moleculesSystem.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <time.h>

using namespace std;

int main(int argc, char **argv) {
	
	srand(time(NULL));
	
//wczytywanie pliku
	parameters param(argv[1]);
	param.readFile();
	
// krawędzie kryształu
	param.setCristalEdges();


//początkowe wartości p i r
	moleculesSystem molSys(argv[2], argv[3] ,param);
	molSys.initMoleculesSystem();
	
//początkowe wartości F i V
	molSys.setMeasures();
	molSys.setSystemMeasures();
	
	molSys.writetHVTPValues();
	
//lab03

	molSys.simulateSystem();
	
	molSys.closeFile();
	
return 0;
}
