#include "JF12.h"

double logisticFunction(double x, double x0, double w) {
	return 1. / (1. + exp(-2. * (fabs(x) - x0) / w));
}

template<typename T>
T getPhi(const T&x, const T&y) {
	T eps = std::numeric_limits<T>::min();
	if ((fabs(x) < eps) and (fabs(y) < eps))
		return 0.0;
	else
		return std::atan2(y, x);
}

JF12Field::JF12Field() {
	useRegular = true;
	useStriated = false;
	useTurbulent = false;

	// spiral arm parameters
	pitch = 11.5 * M_PI / 180;
	sinPitch = std::sin(pitch);
	cosPitch = std::cos(pitch);
	tan90MinusPitch = std::tan(M_PI / 2 - pitch);

	rArms[0] = 5.1 * kpc;
	rArms[1] = 6.3 * kpc;
	rArms[2] = 7.1 * kpc;
	rArms[3] = 8.3 * kpc;
	rArms[4] = 9.8 * kpc;
	rArms[5] = 11.4 * kpc;
	rArms[6] = 12.7 * kpc;
	rArms[7] = 15.5 * kpc;

	// regular field parameters
	bRing = 0.1 * muG;
	hDisk = 0.40 * kpc;
	wDisk = 0.27 * kpc;

	bDisk[0] = 0.1 * muG;
	bDisk[1] = 3.0 * muG;
	bDisk[2] = -0.9 * muG;
	bDisk[3] = -0.8 * muG;
	bDisk[4] = -2.0 * muG;
	bDisk[5] = -4.2 * muG;
	bDisk[6] = 0.0 * muG;
	bDisk[7] = 2.7 * muG;

	bNorth = 1.4 * muG;
	bSouth = -1.1 * muG;
	rNorth = 9.22 * kpc;
	rSouth = 17 * kpc;
	wHalo = 0.20 * kpc;
	z0 = 5.3 * kpc;

	bX = 4.6 * muG;
	thetaX0 = 49.0 * M_PI / 180;
	sinThetaX0 = sin(thetaX0);
	cosThetaX0 = cos(thetaX0);
	tanThetaX0 = tan(thetaX0);
	rXc = 4.8 * kpc;
	rX = 2.9 * kpc;

	// striated field parameter
	sqrtbeta = sqrt(1.36);

	// turbulent field parameters
	bDiskTurb[0] = 10.81 * muG;
	bDiskTurb[1] = 6.96 * muG;
	bDiskTurb[2] = 9.59 * muG;
	bDiskTurb[3] = 6.96 * muG;
	bDiskTurb[4] = 1.96 * muG;
	bDiskTurb[5] = 16.34 * muG;
	bDiskTurb[6] = 37.29 * muG;
	bDiskTurb[7] = 10.35 * muG;

	bDiskTurb5 = 7.63 * muG;
	zDiskTurb = 0.61 * kpc;

	bHaloTurb = 4.68 * muG;
	rHaloTurb = 10.97 * kpc;
	zHaloTurb = 2.84 * kpc;
}

void JF12Field::setUseRegular(bool use) {
	useRegular = use;
}

void JF12Field::setUseStriated(bool use) {
	useStriated = use;
}

void JF12Field::setUseTurbulent(bool use) {
	useTurbulent = use;
}

bool JF12Field::isUsingRegular() {
	return useRegular;
}

bool JF12Field::isUsingStriated() {
	return useStriated;
}

bool JF12Field::isUsingTurbulent() {
	return useTurbulent;
}

std::vector<double> JF12Field::getRegularField(const double& x, const double& y, const double& z) const {
	std::vector<double> b(3, 0);

	double r = sqrt(x * x + y * y); // in-plane radius
	double d = sqrt(x * x + y * y + z * z); // distance to galactic center
	if ((d < 1 * kpc) or (d > 20 * kpc))
		return b; // 0 field for d < 1 kpc or d > 20 kpc

	double phi = getPhi(x, y); // azimuth
	double sinPhi = sin(phi);
	double cosPhi = cos(phi);

	double lfDisk = logisticFunction(z, hDisk, wDisk);

	// disk field
	if (r > 3 * kpc) {
		double bMag;
		if (r < 5 * kpc) {
			// molecular ring
			bMag = bRing * (5 * kpc / r) * (1 - lfDisk);
			b.at(0) += -bMag * sinPhi;
			b.at(1) += bMag * cosPhi;
		} else {
			// spiral region
			double r_negx = r * exp(-(phi - M_PI) / tan90MinusPitch);
			if (r_negx > rArms[7])
				r_negx = r * exp(-(phi + M_PI) / tan90MinusPitch);
			if (r_negx > rArms[7])
				r_negx = r * exp(-(phi + 3 * M_PI) / tan90MinusPitch);

			for (int i = 7; i >= 0; i--)
				if (r_negx < rArms[i])
					bMag = bDisk[i];

			bMag *= (5 * kpc / r) * (1 - lfDisk);
			b.at(0) += bMag * (sinPitch * cosPhi - cosPitch * sinPhi);
			b.at(1) += bMag * (sinPitch * sinPhi + cosPitch * cosPhi);
		}
	}

	// toroidal halo field
	double bMagH = exp(-fabs(z) / z0) * lfDisk;
	if (z >= 0)
		bMagH *= bNorth * (1 - logisticFunction(r, rNorth, wHalo));
	else
		bMagH *= bSouth * (1 - logisticFunction(r, rSouth, wHalo));
	b.at(0) += -bMagH * sinPhi;
	b.at(1) += bMagH * cosPhi;

	// poloidal halo field
	double bMagX;
	double sinThetaX, cosThetaX;
	double rp;
	double rc = rXc + fabs(z) / tanThetaX0;
	if (r < rc) {
		// varying elevation region
		rp = r * rXc / rc;
		bMagX = bX * exp(-1 * rp / rX) * pow(rXc / rc, 2.);
		double thetaX = atan2(fabs(z), (r - rp));
		if (z == 0)
			thetaX = M_PI / 2.;
		sinThetaX = sin(thetaX);
		cosThetaX = cos(thetaX);
	} else {
		// constant elevation region
		rp = r - fabs(z) / tanThetaX0;
		bMagX = bX * exp(-rp / rX) * (rp / r);
		sinThetaX = sinThetaX0;
		cosThetaX = cosThetaX0;
	}
	double zsign = z < 0 ? -1 : 1;
	b.at(0) += zsign * bMagX * cosThetaX * cosPhi;
	b.at(1) += zsign * bMagX * cosThetaX * sinPhi;
	b.at(2) += bMagX * sinThetaX;

	return b;
}

double JF12Field::getTurbulentStrength(const double& x, const double& y, const double& z) const {
	double r = sqrt(x * x + y * y); // in-plane radius
	if (r > 20 * kpc)
		return 0;

	double phi = getPhi(x, y); // azimuth

	// disk
	double bDisk = 0;
	if (r < 5 * kpc) {
		bDisk = bDiskTurb5;
	} else {
		// spiral region
		double r_negx = r * exp(-(phi - M_PI) / tan90MinusPitch);
		if (r_negx > rArms[7])
			r_negx = r * exp(-(phi + M_PI) / tan90MinusPitch);
		if (r_negx > rArms[7])
			r_negx = r * exp(-(phi + 3 * M_PI) / tan90MinusPitch);

		for (int i = 7; i >= 0; i--)
			if (r_negx < rArms[i])
				bDisk = bDiskTurb[i];

		bDisk *= (5 * kpc) / r;
	}
	bDisk *= exp(-0.5 * pow(z / zDiskTurb, 2));

	// halo
	double bHalo = bHaloTurb * exp(-r / rHaloTurb) * exp(-0.5 * pow(z / zHaloTurb, 2));

	// modulate turbulent field
	return sqrt(pow(bDisk, 2) + pow(bHalo, 2));
}

std::vector<double> JF12Field::getTurbulentField(const double& x, const double& y, const double& z) const {
	std::vector<double> b(3, 0); // TODO mettere un campo random
	return b; // (b * getTurbulentStrength(x, y, z));
}

std::vector<double> JF12Field::GetB(const double& x, const double& y, const double& z) {
	std::vector<double> b(3, 0);
	if (useTurbulent) {
		std::vector<double> b_t = getTurbulentField(x, y, z);
		b.at(0) += b_t.at(0);
		b.at(1) += b_t.at(1);
		b.at(2) += b_t.at(2);
	}
	/*if (useStriated) {
	 std::vector<double> b_s = getStriatedField(x, y, z);
	 b.at(0) += b_s.at(0);
	 b.at(1) += b_s.at(1);
	 b.at(2) += b_s.at(2);
	 }*/
	if (useRegular) {
		std::vector<double> b_r = getRegularField(x, y, z);
		b.at(0) += b_r.at(0);
		b.at(1) += b_r.at(1);
		b.at(2) += b_r.at(2);
	}
	return b;
}
