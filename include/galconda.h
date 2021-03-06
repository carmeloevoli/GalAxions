#ifndef __GALCONDA_H
#define __GALCONDA_H

#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <memory>
#include <fstream>
#include <string>
#include <complex>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gas.h"
#include "magneticfield.h"
#include "los.h"
#include "Pshirkov.h"
#include "JF12.h"
#include "mymatrix.h"
#include "solvers.h"

enum MagneticFieldType {
  CONSTANT, FARRAR, PSHIRKOV
};

enum GasDensityType {
  FERRIERE, CORDES, YMW
};

class GALCONDA {
private:

    double axion_mass;
    double g_ag;

    std::string output_filename;
    std::ofstream output_ss;

    std::string los_filename;
    std::ofstream los_ss;

    std::shared_ptr<MagneticField> magnetic_field;
    std::shared_ptr<ElectronModel> gas;

    LOS los;

public:

    GALCONDA(const double& axion_mass_, const double& g_ag_, const std::string& filename_)
            :axion_mass(axion_mass_), g_ag(g_ag_) {
        output_filename = filename_ + ".txt";
        los_filename = filename_ + ".los";
        output_ss.open(output_filename.c_str(), std::ofstream::out);
        los_ss.open(los_filename.c_str(), std::ofstream::out);
    }

    ~GALCONDA() {
        output_ss.close();
        los_ss.close();
    }

    void createGasDensity(const GasDensityType& gastype_);

    void createMagneticField(const MagneticFieldType& btype_, const int& bmode_ = ASS);

    void createLos(const double& ldeg_, const double& bdeg_, const double& maxDistance_ = 30 * kpc);

    void printLos(const double&);

    void calculateProbability(const size_t& nEnergy_, const double& Emin_, const double& Emax_, const bool& do_damping_,
            const bool& do_output_ = true);
};

#endif
