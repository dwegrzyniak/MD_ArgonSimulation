#include "molecule.hpp"
#include <iostream>

void molecule::setInitValues(int i0, int i1, int i2){
	
	coordinates lambda;
	int signx, signy, signz; //zmienne pomocnicze
	
	//std:: cout << (i0 - ((double)param.n-1)/2) * param.b[0].x << '\t' << (i1 - (param.n-1)/2) * param.b[1].x << '\t' << (i2 - (param.n-1)/2) * param.b[2].x << '\t' << std::endl;
	
	r.x = (i0 - ((double)param.n-1)/2) * param.b[0].x + (i1 - (double)(param.n-1)/2) * param.b[1].x + (i2 - (double)(param.n-1)/2) * param.b[2].x;
	r.y = (i0 - ((double)param.n-1)/2) * param.b[0].y + (i1 - (double)(param.n-1)/2) * param.b[1].y + (i2 - (double)(param.n-1)/2) * param.b[2].y;
	r.z = (i0 - ((double)param.n-1)/2) * param.b[0].z + (i1 - (double)(param.n-1)/2) * param.b[1].z + (i2 - (double)(param.n-1)/2) * param.b[2].z;
				
	lambda.x = (double) rand()/RAND_MAX;
	lambda.y = (double) rand()/RAND_MAX;
	lambda.z = (double) rand()/RAND_MAX;
				
	E.x = param.E0 * log(lambda.x);
	E.y = param.E0 * log(lambda.y);
	E.z = param.E0 * log(lambda.z);

	signx = -1;
	signy = -1;
	signz = -1;
		
	if (rand()%2 == 0) signx = 1; 
	if (rand()%2 == 0) signy = 1; 
	if (rand()%2 == 0) signz = 1; 
				
	p.x = (double)signx * pow(2* param.m * E.x, 0.5);
	p.y = (double)signy * pow(2* param.m * E.y, 0.5);
	p.z = (double)signz * pow(2* param.m * E.z, 0.5);
	
}

void molecule::setR() {
	r.x = r.x + (p_tmp.x * param.tau / param.m);
	r.y = r.y + (p_tmp.y * param.tau / param.m);
	r.z = r.z + (p_tmp.z * param.tau / param.m);
	
}

void molecule::setPTmp () {
	p_tmp.x =  p.x + 0.5 * param.tau * (Fs.x + Fij.x);
	p_tmp.y =  p.y + 0.5 * param.tau * (Fs.y + Fij.y);
	p_tmp.z =  p.z + 0.5 * param.tau * (Fs.z + Fij.z);	


}

void molecule::setP() {
	p.x = p_tmp.x + 0.5 * (Fs.x + Fij.x) * param.tau;
	p.y = p_tmp.y + 0.5 * (Fs.y + Fij.y) * param.tau;
	p.z = p_tmp.z + 0.5 * (Fs.z + Fij.z) * param.tau; 
}

molecule::molecule(parameters param0){
	param = param0;
	
}

molecule::molecule(){

}
