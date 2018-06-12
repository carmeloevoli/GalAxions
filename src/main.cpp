#include <iostream>
#include <string>

#include "galconda.h"
#include "cgs.h"

int main(int argc, char** argv) {
    std::cout << "# Welcome to GalAxions!" << std::endl;

    if (argc != 6) {
        std::cerr << "Usage: ./los <output file name> mass [eV] coupling [GeV^-1] l los [deg] b los [deg]\n";
        std::cerr << "Example: ./los HESSJ1640-465 1e-9 1e-10 338.32 -0.02\n";
        return -1;
    }

    const double axion_mass = atof(argv[2]) * eV;
    const double g_ag = atof(argv[3]) / GeV;
    const double l_los_deg = atof(argv[4]);
    const double b_los_deg = atof(argv[5]);

    const std::string filename = argv[1];

    GALCONDA* g = new GALCONDA(axion_mass, g_ag, filename);

    g->createMagneticField(FARRAR); // Alternative models: PSHIRKOV, CONSTANT
    g->createGasDensity(YMW); // Alternative models: CORDES91, FERRIERE
    g->createLos(l_los_deg, b_los_deg);
    g->printLos(30 * kpc);

    g->calculateProbability(1000, 1e-5 * TeV, 1e0 * TeV, false); // HESS J1640-465

    delete g;
    return 0;
}
