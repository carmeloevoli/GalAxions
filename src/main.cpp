#include <iostream>
#include <string>

#include "mks.h"
#include "galaxions.h"

int main(int argc, char** argv) {
	std::cout << "#Welcome to GalAxions!" << std::endl;

	if (argc != 6) {
		std::cerr << "Usage: ./los <output file name> mass [eV] coupling [GeV^-1] l los [deg] b los [deg]\n";
		return -1;
	}

	const double axionMass = atof(argv[2]) * eV;
	const double gag = atof(argv[3]) / GeV;
	const double l_los = atof(argv[4]);
	const double b_los = atof(argv[5]);

	const std::string initFilename = argv[1];

	galAxions * g = new galAxions(axionMass, gag, initFilename);

	g->createMagneticField(FARRAR); // Alternative models: PSHIRKOV, CONSTANT
	g->createGasDensity(FERRIERE);
	g->createLos(l_los, b_los);
	g->printLos(10 * kpc);

	//g->calculateProbability(100,1.0*KEV,1e2*KEV,false,true);
	//g->calculateProbability(2000,1.0*MEV,100.0*MEV,false,true); // Payez Model
	//g->calculateProbability(2000,1.0*MEV,300.0*MEV,false,true); // Payez Model GC
	//g->calculateProbability(2000, 1e-4 * TEV, 1e2 * TEV, false, true); // HESS J1640-465

	if (g)
		delete g;
	return 0;
}
