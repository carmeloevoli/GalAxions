#include "Ferriere07.h"

double phi(const double& r_pc, const double& z_pc) {
	double C1 = 8.887e3;
	double C2 = 3.e3;
	double C3 = 0.325;
	double a1 = 6.5e3;
	double a2 = 0.7e3;
	double a3 = 12e3;
	double b1 = 0.26e3;
	double rh = 210e3;
	double P1 = C1 / std::sqrt(pow2(r_pc) + pow2(a1 + std::sqrt(pow2(z_pc) + pow2(b1))));
	double P2 = C2 / (a2 + std::sqrt(pow2(r_pc) + pow2(z_pc)));
	double L = std::sqrt(1. + (pow2(a3) + pow2(r_pc) + pow2(z_pc)) / pow2(rh));
	double P3 = -C3 * std::log((L - 1) / (L + 1));
	return -pow2(225.) * (P1 + P2 + P3);
}

double Ferriere07::wim(double x_pc, double y_pc, double z_pc) const {
	double r = std::sqrt(pow2(x_pc) + pow2(y_pc));
	double y_3 = -10;
	double z_3 = -20;

	double P1 = exp(-(pow2(x_pc) + pow2(y_pc - y_3)) / pow2(L3)) * exp(-(pow2(z_pc - z_3)) / pow2(H3));

	double P2 = 0.009;
	P2 *= exp(-(pow2(r - L2)) / (pow2(L2) / 4.));
	P2 *= pow2(sech(z_pc / H2));

	double P3 = 0.005;
	P3 *= std::cos(M_PI * r / 2. / L1);
	P3 *= pow2(sech(z_pc / H1));

	return 8.0 * (P1 + P2 + P3);
}

double Ferriere07::vhim(double x_pc, double y_pc, double z_pc) const {
	double eta = y_pc * cos_alphavh + z_pc * sin_alphavh;
	double zeta = -y_pc * sin_alphavh + z_pc * cos_alphavh;
	double A = (pow2(x_pc) + pow2(eta)) / pow2(Lvh);
	double B = pow2(zeta) / pow2(Hvh);
	return 0.29 * std::exp(-(A + B));
}

double Ferriere07::him(double x_pc, double y_pc, double z_pc) const {
	double r_pc = std::sqrt(x_pc * x_pc + y_pc * y_pc);
	double km2cm = 1e5;
	double delta_phi = (phi(r_pc, z_pc) - phi(0., 0.)) * pow2(km2cm);
	return std::pow(std::pow(0.009, 2. / 3.) - 1.54e-17 * delta_phi, 1.5);
}

double Ferriere07::disk(const double& r, const double& z) const {
	double value = 2.13e-2 / cm3;
	value *= std::exp(-(r - sun_r) / (18.24 * kpc));
	value *= std::min(1.0, std::exp((fabs(z) - 1.0 * kpc) / (0.5 * kpc)));
	return value;
}

double Ferriere07::bulge(const double& x, const double& y, const double& z) const {
	double value = 0;
	if (fabs(x) < kpc && fabs(y) < kpc && fabs(z) < kpc) {
		value += wim(x / pc, y / pc, z / pc);
		value += vhim(x / pc, y / pc, z / pc);
		value += him(x / pc, y / pc, z / pc);
	}
	return value / cm3;
}
