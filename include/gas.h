#ifndef GALCONDA_GAS_DISTRIBUTION_H_
#define GALCONDA_GAS_DISTRIBUTION_H_

#include <algorithm>
#include <cmath>

#include "cgs.h"
#include "constants.h"

extern "C" double ymw16_ne(double x_pc, double y_pc, double z_pc);

class ElectronModel {
public:
    ElectronModel() {
    }

    virtual ~ElectronModel() {
    }

    virtual double get(const double& x, const double& y, const double& z) const = 0;
};

class Cordes91 : public ElectronModel {
private:
    double fne1 = 0.025 / cm3;
    double H1 = 1.00 * kpc;
    double A1 = 20.0 * kpc;
    double fne2 = 0.200 / cm3;
    double H2 = 0.15 * kpc;
    double A2 = 2.0 * kpc;
    double R2 = 4.0 * kpc;

public:
    double get(const double& x, const double& y, const double& z) const override;
};

class YMW16 : public ElectronModel {
public:
    double get(const double& x, const double& y, const double& z) const override;
};

#endif
