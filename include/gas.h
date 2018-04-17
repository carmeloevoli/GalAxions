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
