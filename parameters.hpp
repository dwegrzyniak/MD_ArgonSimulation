#ifndef _parameters_hpp_
#define _parameters_hpp_

#include <fstream>
#include <string>
#include <cmath>

struct coordinates {
	double x;
	double y;
	double z;
};

class parameters{
public:	
	std::string fileName;
	
	double k = 8.31e-3;
	double pi = 3.14159;
	
	int n, N;
	double m, epsilon, tau, f, a, R, L, T0, E0; 
	int S_o, S_d, S_out, S_xyz;
	
	coordinates b[3];

public:
	parameters();
	parameters(std::string file);
	void readFile();
	void setCristalEdges();

};
#endif
