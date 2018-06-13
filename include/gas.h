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

class Ferriere07 : public ElectronModel {
private:
    double L1 = 17000;
    double H1 = 950;
    double L2 = 3700;
    double H2 = 140;
    double L3 = 145;
    double H3 = 26;
    double Lvh = 162;
    double Hvh = 90;
    double alpha_vh = deg2rad(21.);
    double cos_alpha_vh = std::cos(alpha_vh);
    double sin_alpha_vh = std::sin(alpha_vh);

    double wim(double x_pc, double y_pc, double z_pc) const;

    double vhim(double x_pc, double y_pc, double z_pc) const;

    double him(double x_pc, double y_pc, double z_pc) const;

public:
    double disk(const double& r, const double& z) const;

    double bulge(const double& x, const double& y, const double& z) const;

    double get(const double& x, const double& y, const double& z) const override {
        double r = std::sqrt(x * x + y * y);
        return std::max(disk(r, z) + bulge(x, y, z), 0.);
    }
};

#endif
