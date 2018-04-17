#include "gas.h"

double phi(const double& r, const double& z) {
	double C1 = 8.887e3;
	double C2 = 3.e3;
	double C3 = 0.325;
	double a1 = 6.5e3;
	double a2 = 0.7e3;
	double a3 = 12e3;
	double b1 = 0.26e3;
	double rh = 210e3;
	double P1 = C1 / std::sqrt(pow2(r) + pow2(a1 + std::sqrt(pow2(z) + pow2(b1))));
	double P2 = C2 / (a2 + std::sqrt(pow2(r) + pow2(z)));
	double L = std::sqrt(1. + (pow2(a3) + pow2(r) + pow2(z)) / pow2(rh));
	double P3 = -C3 * std::log((L - 1) / (L + 1));
	return -pow2(225.) * (P1 + P2 + P3);
}

double Ferriere07::wim(double x, double y, double z) const {
	double r = std::sqrt(pow2(x) + pow2(y));
	double y_3 = -10;
	double z_3 = -20;

	double P1 = exp(-(pow2(x) + pow2(y - y_3)) / pow2(L3)) * exp(-(pow2(z - z_3)) / pow2(H3));

	double P2 = 0.009;
	P2 *= exp(-(pow2(r - L2)) / (pow2(L2) / 4.));
	P2 *= pow2(sech(z / H2));

	double P3 = 0.005;
	P3 *= std::cos(M_PI * r / 2. / L1);
	P3 *= pow2(sech(z / H1));

	return 8.0 * (P1 + P2 + P3);
}

double Ferriere07::vhim(double x, double y, double z) const {
	double eta = y * cos_alphavh + z * sin_alphavh;
	double zeta = -y * sin_alphavh + z * cos_alphavh;
	double A = (pow2(x) + pow2(eta)) / pow2(Lvh);
	double B = pow2(zeta) / pow2(Hvh);
	return 0.29 * std::exp(-(A + B));
}

double Ferriere07::him(double x, double y, double z) const {
	double r = std::sqrt(x * x + y * y);
	double km2cm = 1e5;
	double delta_phi = (phi(r, z) - phi(0., 0.)) * pow2(km2cm);
	return std::pow(std::pow(0.009, 2. / 3.) - 1.54e-17 * delta_phi, 1.5);
}

