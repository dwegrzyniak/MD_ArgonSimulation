#ifndef _molecule_hpp_
#define _molecule_hpp_

#include "parameters.hpp"

#include <time.h>

class molecule {
	
public:
	parameters param;

	coordinates r;
	coordinates E;
	coordinates p, p_tmp;
	
	double Vs, Vp; //potencjały
	coordinates Fs, Fij;
public:
	molecule();
	molecule(parameters param0);
	void setInitValues(int i0, int i1, int i2);
	void setPTmp();
	void setR();
	void setP();
};



#endif
