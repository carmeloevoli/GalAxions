#ifndef __MAGNFIELD_H
#define __MAGNFIELD_H

#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>

#include "mks.h"
#include "constants.h"

inline double normVector(const std::vector<double>& v) {
	return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

class MagneticField {

private:
	double Bmin = 0;
	double Bmax = 0;

public:
	MagneticField() {
	}

	virtual ~MagneticField() {
	}

	std::vector<double> GetBperp(const double&, const double&, const double&);

	double GetBtransverse(const double&, const double&, const double&);

	virtual std::vector<double> GetB(const double&, const double&, const double&) {
		return std::vector<double>(3, 0);
	}
};

class ConstantField: public MagneticField {

private:
	double BConstant = 1. * muG;

public:
	ConstantField() {
	}

	ConstantField(const double& B_) {
		BConstant = B_;
	}

	virtual ~ConstantField() {
	}

	virtual std::vector<double> GetB(const double& x, const double& y, const double& z) {
			std::vector<double> Bret;
			Bret.push_back(BConstant / std::sqrt(3.0));
			Bret.push_back(BConstant / std::sqrt(3.0));
			Bret.push_back(BConstant / std::sqrt(3.0));
			return Bret;
	}
};

#endif
