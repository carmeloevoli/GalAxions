#ifndef _PSHIRKOVFIELD_H
#define _PSHIRKOVFIELD_H

#include <cmath>
#include <vector>
#include <cstdlib>
#include "magneticfield.h"
#include "constants.h"

#define ASS 0
#define BSS 1

class PshirkovField: public MagneticField {

private:
	unsigned int bMode = 0; // [ASS,BSS]
	double p = 0; // [ASS,BSS]
	double z0 = 0;
	double d = 0;
	double B0 = 0;
	double Rc = 0;
	double z0H = 0;
	double R0H = 0;
	double B0H[2]; // [mug] [ASS south, ASS north+BSS]
	double z1H[2]; // [kpc] [inner,outer]

public:
	PshirkovField() {
	}
	PshirkovField(const unsigned int&, const double&, const double&, const double&);
	virtual ~PshirkovField() {
	}

	virtual std::vector<double> GetB(const double&, const double&, const double&);

	double GetBdisk(const double&, const double&, const double&);
	double GetBhalo(const double&, const double&);
};

#endif
