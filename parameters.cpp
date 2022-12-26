#include "parameters.hpp"

using namespace std;

parameters::parameters(string file){
	fileName = file;	
}

parameters::parameters(){
	
}

void parameters::readFile(){
	ifstream parametersFile;
	parametersFile.open(fileName);
	
	parametersFile >> n >> m >> epsilon >> R >> f >> L >> a >> T0 >> tau >> S_o >> S_d  >> S_out >> S_xyz;
	N = n*n*n;
	E0 = -0.5 * T0 * k;
	
	parametersFile.close();
}

void parameters::setCristalEdges(){
	b[0].x = a;
	b[0].y = 0;
	b[0].z = 0;
	
	b[1].x = a/2;
	b[1].y = a * pow(3, 0.5)/2;
	b[1].z = 0;
	
	b[2].x = a/2;
	b[2].y = a * pow(3, 0.5)/6;
	b[2].z = a * pow((double)2/3, 0.5);
}
