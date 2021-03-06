#ifndef __SOLVERS_H
#define __SOLVERS_H

#include <cstdlib>
#include <vector>
#include <iostream>

#include "los.h"
#include "mymatrix.h"
#include "constants.h"

void FairbairnSolver(const int& Ndom, const std::vector<double>& DeltaAgamma, const std::vector<double>& DeltaPl,
        const double& DeltaA, const std::vector<double>& psik, const std::vector<double>& gas_dens,
        const std::vector<double>&, const double&, const double&, const MyMatrix& InitialCondition, MyMatrix& Rho);

void Solver(const int& Ndom, const std::vector<double>& DeltaAgamma, const std::vector<double>& DeltaPl,
        const std::vector<double>& DeltaQED, const std::vector<double>& DeltaPar, const std::vector<double>& DeltaPerp,
        const double& DeltaA, const std::vector<domain>& los, const double& xsecH2, const double& xsecHI,
        const MyMatrix& OldRho, MyMatrix& Rho, MyMatrix& Tktotal);

#endif
