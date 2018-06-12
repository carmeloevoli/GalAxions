#include "solvers.h"

void Solver(const int& Ndom, const std::vector<double>& DeltaAgamma, const std::vector<double>& DeltaPl,
        const std::vector<double>& DeltaQED, const std::vector<double>& DeltaPar, const std::vector<double>& DeltaPerp,
        const double& DeltaA, const std::vector<domain>& domains, const double& xsecH2, const double& xsecHI,
        const MyMatrix& OldRho, MyMatrix& Rho, MyMatrix& Tktotal) {

    const double L = step_size; // [kpc]

    MyMatrix TA, TB, TC;
    MyMatrix add1, add2, add3;
    MyMatrix Tkcross, Tktemp, Tkfinal(true), Tkfinaltemp;
    MyMatrix Tkone(true);

    double D[3];

    for (int i = 0; i < Ndom; ++i) {

        const double cp = std::cos(domains.at(i).psik);
        const double sp = std::sin(domains.at(i).psik);

        const double Gamma =
                3.0 * (domains.at(i).H2_density * xsecH2 + domains.at(i).H2_density * xsecHI); // TODO CHECK!

        const double Theta = 0.5 * atan2(2.0 * DeltaAgamma[i], (DeltaPar[i] - DeltaA)); // Eq. 3.19

        InitD(DeltaPerp[i], DeltaPar[i], DeltaA, DeltaAgamma[i], D); // Eq. 3.16

        InitTABC(cp, sp, Theta, TA, TB, TC); // Eq. 3.38

        Mult(TA, exp(std::complex<double>(-0.5 * Gamma * L, D[0] * L)), add1);
        Mult(TB, exp(std::complex<double>(-0.5 * Gamma * L, D[1] * L)), add2);
        Mult(TC, exp(std::complex<double>(0.0, D[2] * L)), add3);

        Add(add1, add2, add3, Tktemp);

        Mult(Tkone, Tktemp, Tktotal);
        Conjugate(Tktotal, Tkcross);

        Mult(Tktotal, OldRho, Tkcross, Rho);
        Equal(Tktotal, Tkone);

#ifdef DEBUGMODE
        auto Pag_PDF_ = real(Rho(2,2));
        auto Iav_PDF_ = real(Rho(0,0) + Rho(1,1));
        std::cout << Iav_PDF_ << "\t" << Pag_PDF_ << "\n";
#endif
    }
}
