#include "solvers.h"

void BassanSolver(const int& Ndom, const std::vector<double>& DeltaAgamma, const std::vector<double>& DeltaPl, const std::vector<double>& DeltaQED,
		const std::vector<double>& DeltaPar, const std::vector<double>& DeltaPerp, const double& DeltaA, const std::vector<double>& psik,
		const std::vector<double>& gas_densH2, const std::vector<double>& gas_densHI, const double& xsecH2, const double& xsecHI, const MyMatrix& OldRho,
		MyMatrix& Rho, MyMatrix& Tktotal) {

	const double L = step_size; // [kpc]

	MyMatrix TA, TB, TC;
	MyMatrix add1, add2, add3;
	MyMatrix Tkcross, Tktemp, Tkfinal(true), Tkfinaltemp;
	MyMatrix Tkone(true);

	double cp, sp;

	double D[3];

	for (int ione = 0; ione < Ndom; ++ione) {

		cp = cos(psik[ione]);
		sp = sin(psik[ione]);

		const double Gamma = 3.0 * (gas_densH2[ione] * xsecH2 + gas_densHI[ione] * xsecHI); // CHECK!

		const double Theta = 0.5 * atan2(2.0 * DeltaAgamma[ione], (DeltaPar[ione] - DeltaA)); // Eq. 3.19

		InitD(DeltaPerp[ione], DeltaPar[ione], DeltaA, DeltaAgamma[ione], D); // Eq. 3.16

		InitTABC(cp, sp, Theta, TA, TB, TC); // Eq. 3.38

		Mult(TA, exp(std::complex<double>(-0.5 * Gamma * L, D[0] * L)), add1);
		Mult(TB, exp(std::complex<double>(-0.5 * Gamma * L, D[1] * L)), add2);
		Mult(TC, exp(std::complex<double>(0.0, D[2] * L)), add3);

		Add(add1, add2, add3, Tktemp);

		Mult(Tkone, Tktemp, Tktotal);
		Conjugate(Tktotal, Tkcross);

		/*
		 if (PRINT_AT_EACH_TIMESTEP) {
		 Mult(Tktotal,Tkfinal,Tkfinaltemp);
		 Equal(Tkfinaltemp,Tkfinal);
		 }
		 */

		Mult(Tktotal, OldRho, Tkcross, Rho);
		Equal(Tktotal, Tkone);

#ifdef DEBUGMODE
		const double PagPDF_temp = real(Rho(2,2));
		const double IavPDF_temp = real(Rho(0,0) + Rho(1,1));
		std::cout<<IavPDF_temp<<" "<<PagPDF_temp<<std::endl;
#endif
	}

	return;
}
