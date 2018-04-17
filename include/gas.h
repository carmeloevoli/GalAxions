#ifndef __GAS_DISTRIBUTION_H
#define __GAS_DISTRIBUTION_H

#include <algorithm>
#include <cmath>
#include "constants.h"

class ElectronModel {
public:
	ElectronModel() {
	}
	virtual ~ElectronModel() {
	}
	virtual double get(const double& x, const double& y, const double& z) const = 0;
};

class Ferriere07: public ElectronModel {
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
	double wim(double x, double y, double z) const;
	double vhim(double x, double y, double z) const;
	double him(double x, double y, double z) const;
public:
	double disk(const double& r, const double& z) const {
		double value = 2.13e-2 / cm3;
		value *= std::exp(-(r - sun_r) / (18.24 * kpc));
		value *= std::min(1.0, std::exp((fabs(z) - 1.0 * kpc) / (0.5 * kpc)));
		return value;
	}

	double bulge(const double& x, const double& y, const double& z) const {
		double value = 0;
		value += wim(x / pc, y / pc, z / pc);
		value += vhim(x / pc, y / pc, z / pc);
		value += him(x / pc, y / pc, z / pc);
		return value / cm3;
	}

	double get(const double& x, const double& y, const double& z) const override {
		double r = std::sqrt(x * x + y * y);
		return std::max(disk(r, z) + bulge(x, y, z), 0.);
	}
};

class Cordes91: public ElectronModel {
private:
	double fne1 = 0.025 / cm3;
	double H1 = 1.00 * kpc;
	double A1 = 20.0 * kpc;
	double fne2 = 0.200 / cm3;
	double H2 = 0.15 * kpc;
	double A2 = 2.0 * kpc;
	double R2 = 4.0 * kpc;
public:
	double get(const double& x, const double& y, const double& z) const override {
		double r = std::sqrt(x * x + y * y);
		double ne1 = fne1 * exp(-fabs(z) / H1) * exp(-pow2(r / A1));
		double ne2 = fne2 * exp(-fabs(z) / H2) * exp(-pow2((r - R2) / A2));
		return ne1 + ne2;
	}
};

class YMW16: public ElectronModel {
public:
	double get(const double& x, const double& y, const double& z) const override {
		//return ymw16_ne(pos.y / pc, pos.x / pc, pos.z / pc) / cm3;
		return 0;
	}
};

#endif
