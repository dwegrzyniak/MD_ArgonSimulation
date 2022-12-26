#ifndef _moleculesSystem_hpp_
#define _moleculesSystem_hpp_

#include "parameters.hpp"
#include "molecule.hpp"
#include "parameters.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <cmath>

class moleculesSystem {
public:
	std::string xyzFileName, outFileName;
	std::ofstream xyzFile, outFile;
	std::vector <molecule> molecules;
	parameters param;
	double V, P, t, H, T;
	
public:
	moleculesSystem(std::string xyzfile, std::string outfile, parameters param0);
	void initMoleculesSystem();
	void setMeasures();
	void simulateSystem();
	void setSystemMeasures();
	void writetHVTPValues(); //  H, V, T, P
	void writexyzeValues();
	void resetForces();
	void resetVp();
	void closeFile();
};


#endif
