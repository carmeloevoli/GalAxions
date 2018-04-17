#include "magneticfield.h"

std::vector<double> MagneticField::GetBperp(const double& x, const double& y, const double& z) {

	const double dist = std::sqrt(pow2(x - sun_x) + pow2(y - sun_y) + pow2(z - sun_z));

	std::vector<double> Bperp(3, 0.0);
	std::vector<double> Btotal = GetB(x, y, z);
	std::vector<double> Versor(3, 0.0);

	Versor[0] = (x - sun_x) / dist;
	Versor[1] = (y - sun_y) / dist;
	Versor[2] = (z - sun_z) / dist;

	Bperp[0] = -(Btotal[1] * Versor[2] - Btotal[2] * Versor[1]);
	Bperp[1] = (Btotal[0] * Versor[2] - Btotal[2] * Versor[0]);
	Bperp[2] = -(Btotal[0] * Versor[1] - Btotal[1] * Versor[0]);

#ifdef DEBUGMODE
	double BtotalDotVersor = Btotal[0] * Versor[0] + Btotal[1] * Versor[1] + Btotal[2] * Versor[2];
	double BperpDotVersor = Bperp[0] * Versor[0] + Bperp[1] * Versor[1] + Bperp[2] * Versor[2];
	double BtotalNorm = normVector(Btotal);
	double BperpNorm = normVector(Bperp);
	double VersorNorm = normVector(Versor);
	std::cout << std::scientific;
	std::cout << Versor[0] << "\t" << Versor[1] << "\t" << Versor[2] << "\t" << Bperp[0] << "\t" << Bperp[1] << "\t" << Bperp[2] << "\t";
	std::cout << BtotalDotVersor / BtotalNorm << "\t" << std::sin(std::acos(BtotalDotVersor / BtotalNorm)) << "\t" << BperpNorm / BtotalNorm << std::endl;
#endif

	return Bperp;
}

