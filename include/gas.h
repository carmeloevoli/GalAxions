#ifndef __GAS_DISTRIBUTION_H
#define __GAS_DISTRIBUTION_H

#include <cmath>
#include "constants.h"

class ElectronModel {
public:
	ElectronModel() {
	}
	virtual ~ElectronModel() {
	}
	virtual double get(const double& r, const double& z) const = 0;
	virtual double get(const double& x, const double& y, const double& z) const = 0;
};

class Ferriere : public ElectronModel {
public:
	double get(const double& r, const double& z) const override {
		double value = 2.13e-2 / cm3;
		value *= std::exp(-(r - sun_r) / (18.24 * kpc));
		value *= std::min(1.0, std::exp((fabs(z) - 1.0 * kpc) / (0.5 * kpc)));
		return value;
	}

	double get(const double& x, const double& y, const double& z) const override {
		return 0;
	}
};

#endif
