#include "moleculesSystem.hpp"

using namespace std;

moleculesSystem::moleculesSystem(string xyzfile, string outfile, parameters param0){
	param = param0;
	xyzFile.open(xyzfile);
	outFile.open(outfile);
	t = 0;
	T = param.T0;
	molecules.resize(param.N);
}

void moleculesSystem::simulateSystem(){
	
	for (int s = 1; s <= param.S_o + param.S_d; s++) {
		
		t += param.tau;
		
		for (int i = 0; i < param.N; i++) {
			molecules[i].setPTmp();
			molecules[i].setR();
		}
		
		setMeasures();
		
		for (int i = 0; i < param.N; i++) {
			molecules[i].setP();
		}
		
		setSystemMeasures();
		
		if (s % param.S_out == 0 || t == 0) {
			writetHVTPValues();
		}
		
		if (s % param.S_xyz == 0) {
			writexyzeValues();
		}
		
	}
	
}

void moleculesSystem::setSystemMeasures(){
	H = 0;
	T = 0;
	V = 0;
	P = 0;
	for (int i = 0; i < param.N; i++) {
		double p_square = molecules[i].p.x * molecules[i].p.x + molecules[i].p.y * molecules[i].p.y + molecules[i].p.z * molecules[i].p.z;
		
		H +=p_square /(2 * param.m);
		
		T +=(2 * p_square)/(2 * 3 * param.N * param.k * param.m);
		
		V+= molecules[i].Vs;
		V+= molecules[i].Vp;
		
	//(15)
		P += pow(molecules[i].Fs.x * molecules[i].Fs.x + molecules[i].Fs.y * molecules[i].Fs.y + molecules[i].Fs.z * molecules[i].Fs.z, 0.5)/ (4 * param.pi * param.L * param.L);
	}
	
	H += V;
}

void moleculesSystem::resetForces () {
	for (int i = 0; i < param.N; i++){
		molecules[i].Fij.x = 0;
		molecules[i].Fij.y = 0;
		molecules[i].Fij.z = 0;
	}
	
}

void moleculesSystem::resetVp () {
	for (int i = 0; i < param.N; i++){
		molecules[i].Vp = 0;
	}
	
}

void moleculesSystem::setMeasures(){
	double ri, rij;
	resetVp();
	resetForces();
	for (int i = 0; i < param.N; i++) {
	//(10)
		ri = pow(molecules[i].r.x * molecules[i].r.x + molecules[i].r.y * molecules[i].r.y + molecules[i].r.z * molecules[i].r.z, 0.5);

		if (ri < param.L) molecules[i].Vs = 0;
		else
		{
			molecules[i].Vs = 0.5 * param.f * (ri - param.L) * (ri - param.L);
		}
		
	//(14)
		if (ri < param.L) {
			 molecules[i].Fs.x = 0;
			 molecules[i].Fs.y = 0;
			 molecules[i].Fs.z = 0;
		 }
		else {
			molecules[i].Fs.x = param.f * (param.L - ri) * molecules[i].r.x / ri;
			molecules[i].Fs.y = param.f * (param.L - ri) * molecules[i].r.y / ri;
			molecules[i].Fs.z = param.f * (param.L - ri) * molecules[i].r.z / ri;
		}
		
		
		for (int j = 0; j < i; j++){
		//(9)
			rij = pow( pow( molecules[i].r.x - molecules[j].r.x, 2) + pow( molecules[i].r.y - molecules[j].r.y, 2) + pow( molecules[i].r.z - molecules[j].r.z, 2), 0.5);
			molecules[i].Vp += param.epsilon * ( pow(param.R/rij, 12) - 2 * pow(param.R/rij, 6) );
			
		//(13)
			molecules[i].Fij.x += 12 * param.epsilon * ( pow(param.R/rij, 12) - pow(param.R/rij, 6)) * (molecules[i].r.x - molecules[j].r.x) / (rij * rij);
			molecules[i].Fij.y += 12 * param.epsilon * ( pow(param.R/rij, 12) - pow(param.R/rij, 6)) * (molecules[i].r.y - molecules[j].r.y) / (rij * rij);
			molecules[i].Fij.z += 12 * param.epsilon * ( pow(param.R/rij, 12) - pow(param.R/rij, 6)) * (molecules[i].r.z - molecules[j].r.z) / (rij * rij);
		
			molecules[j].Fij.x -= 12 * param.epsilon * ( pow(param.R/rij, 12) - pow(param.R/rij, 6)) * (molecules[i].r.x - molecules[j].r.x) / (rij * rij);
			molecules[j].Fij.y -= 12 * param.epsilon * ( pow(param.R/rij, 12) - pow(param.R/rij, 6)) * (molecules[i].r.y - molecules[j].r.y) / (rij * rij);
			molecules[j].Fij.z -= 12 * param.epsilon * ( pow(param.R/rij, 12) - pow(param.R/rij, 6)) * (molecules[i].r.z - molecules[j].r.z) / (rij * rij);
		

		}
	}
}

void moleculesSystem::initMoleculesSystem(){
	
	double px = 0;
	double py = 0;
	double pz = 0;
	
	int i = -1;
	
	for (int i0 = 0; i0 < param.n; i0++){
		
		for (int i1 = 0; i1 < param.n; i1++){
			
			for (int i2 = 0; i2 < param.n; i2++){
				i++;
				
				molecule mol(param);
				
				mol.setInitValues(i0, i1, i2);
				molecules[i] = mol;
				px += mol.p.x;
				py += mol.p.y;
				pz += mol.p.z;
				
			}
		}
	}
	
	px = px/param.N;
	py = py/param.N;
	pz = pz/param.N;
	
	
	for (int i = 0; i < param.N; i++){
		molecules[i].p.x -= px;
		molecules[i].p.y -= py;
		molecules[i].p.z -= pz;
	}
	
	
	writexyzeValues();
}

void moleculesSystem::writetHVTPValues(){
	outFile << t << "\t" << H << "\t" << V << "\t" << T << "\t" << P << endl;
}

void moleculesSystem::writexyzeValues(){
	xyzFile << param.N << endl << endl;
	for (int i = 0; i < param.N; i++){
		//double p_square = molecules[i].p.x * molecules[i].p.x + molecules[i].p.y * molecules[i].p.y + molecules[i].p.z * molecules[i].p.z;
		xyzFile << "Ar" << "\t" << molecules[i].r.x << "\t" << molecules[i].r.y << "\t" << molecules[i].r.z << endl;
	}
}

void moleculesSystem::closeFile(){
	xyzFile.close();
	outFile.close();
}
