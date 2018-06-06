#include <iostream>
#include <string>

#include "mks.h"
#include "galaxions.h"

int main(int argc, char** argv) {
	std::cout << "# Welcome to GalAxions!" << std::endl;

	if (argc != 6) {
		std::cerr << "Usage: ./los <output file name> mass [eV] coupling [GeV^-1] l los [deg] b los [deg]\n";
		std::cerr << "Example: ./los HESSJ1640-465 1e-9 1e-10 338.32 -0.02\n";
		return -1;
	}

	const double axionMass = atof(argv[2]) * eV;
	const double gag = atof(argv[3]) / GeV;
	const double l_los_deg = atof(argv[4]);
	const double b_los_deg = atof(argv[5]);

	const std::string initFilename = argv[1];

	galAxions * g = new galAxions(axionMass, gag, initFilename);

	g->createMagneticField(FARRAR); // Alternative models: PSHIRKOV, CONSTANT
	g->createGasDensity(YMW); // Alternative models: CORDES91, FERRIERE
	g->createLos(l_los_deg, b_los_deg);
	g->printLos(30 * kpc);

	g->calculateProbability(200, 1e-4 * TeV, 1e2 * TeV, false); // HESS J1640-465

	if (g)
		delete g;
	return 0;
}
