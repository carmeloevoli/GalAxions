#ifndef GALCONDA_FERRIERE07_H_
#define GALCONDA_FERRIERE07_H_

#include <cmath>
#include "cgs.h"
#include "gas.h"

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
    double alphavh = deg2rad(21.);
    double cos_alphavh = std::cos(alphavh);
    double sin_alphavh = std::sin(alphavh);

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

#endif /* INCLUDE_FERRIERE07_H_ */
