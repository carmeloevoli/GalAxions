#include "mymatrix.h"

void InitD(const double DeltaPerp, const double DeltaPar, const double DeltaA, const double DeltaAg, double D[]) {

    D[0] = DeltaPerp; // Eq. 3.20

    double dsym = 0.5 * (DeltaPar + DeltaA);
    double dasym = 0.5 * sqrt(pow(DeltaA - DeltaPar, 2.0) + 4.0 * pow(DeltaAg, 2.0));

    D[1] = dsym + dasym; // Eq. 3.21
    D[2] = dsym - dasym; // Eq. 3.22
}

void InitTABC(const double cp, const double sp, const double theta, MyMatrix &TA, MyMatrix &TB, MyMatrix &TC) {

    double ct = cos(theta);
    double st = sin(theta);

    double cpcp = cp * cp;
    double stst = st * st;
    double ctct = ct * ct;
    double stct = st * ct;
    double spsp = sp * sp;
    double spcp = sp * cp;

    // Eq. 3.38
    TA(0, 0) = std::complex<double>(cpcp, 0.0);
    TA(0, 1) = std::complex<double>(-spcp, 0.0);
    TA(0, 2) = std::complex<double>(0.0, 0.0);

    TA(1, 0) = TA(0, 1);
    TA(1, 1) = std::complex<double>(spsp, 0.0);
    TA(1, 2) = std::complex<double>(0.0, 0.0);

    TA(2, 0) = TA(0, 2);
    TA(2, 1) = TA(1, 2);
    TA(2, 2) = std::complex<double>(0.0, 0.0);

    // Eq. 3.39
    TB(0, 0) = std::complex<double>(spsp * ctct, 0.0);
    TB(0, 1) = std::complex<double>(spcp * ctct, 0.0);
    TB(0, 2) = std::complex<double>(stct * sp, 0.0);

    TB(1, 0) = TB(0, 1);
    TB(1, 1) = std::complex<double>(cpcp * ctct, 0.0);
    TB(1, 2) = std::complex<double>(stct * cp, 0.0);

    TB(2, 0) = TB(0, 2);
    TB(2, 1) = TB(1, 2);
    TB(2, 2) = std::complex<double>(stst, 0.0);

    // Eq. 3.40
    TC(0, 0) = std::complex<double>(stst * spsp, 0.0);
    TC(0, 1) = std::complex<double>(stst * spcp, 0.0);
    TC(0, 2) = std::complex<double>(-stct * sp, 0.0);

    TC(1, 0) = TC(0, 1);
    TC(1, 1) = std::complex<double>(stst * cpcp, 0.0);
    TC(1, 2) = std::complex<double>(-stct * cp, 0.0);

    TC(2, 0) = TC(0, 2);
    TC(2, 1) = TC(1, 2);
    TC(2, 2) = std::complex<double>(ctct, 0.0);
}

void InitTA(const double cp, const double sp, MyMatrix &TA) {

    TA(0, 0) = std::complex<double>(cp * cp, 0.0);
    TA(0, 1) = std::complex<double>(-sp * cp, 0.0);
    TA(0, 2) = std::complex<double>(0.0, 0.0);

    TA(1, 0) = TA(0, 1);
    TA(1, 1) = std::complex<double>(sp * sp, 0.0);
    TA(1, 2) = std::complex<double>(0.0, 0.0);

    TA(2, 0) = TA(0, 2);
    TA(2, 1) = TA(1, 2);
    TA(2, 2) = std::complex<double>(0.0, 0.0);
}

void InitTB(const double cp, const double sp, const double theta, MyMatrix &TB) {

    double ct = cos(theta);
    double st = sin(theta);

    TB(0, 0) = std::complex<double>(sp * sp * ct * ct, 0.0);
    TB(0, 1) = std::complex<double>(sp * cp * ct * ct, 0.0);
    TB(0, 2) = std::complex<double>(st * ct * sp, 0.0);

    TB(1, 0) = TB(0, 1);
    TB(1, 1) = std::complex<double>(cp * cp * ct * ct, 0.0);
    TB(1, 2) = std::complex<double>(st * ct * cp, 0.0);

    TB(2, 0) = TB(0, 2);
    TB(2, 1) = TB(1, 2);
    TB(2, 2) = std::complex<double>(st * st, 0.0);
}

void InitTC(const double cp, const double sp, const double theta, MyMatrix &TC) {

    double ct = cos(theta);
    double st = sin(theta);

    TC(0, 0) = std::complex<double>(st * st * sp * sp, 0.0);
    TC(0, 1) = std::complex<double>(st * st * sp * cp, 0.0);
    TC(0, 2) = std::complex<double>(-st * ct * sp, 0.0);

    TC(1, 0) = TC(0, 1);
    TC(1, 1) = std::complex<double>(st * st * cp * cp, 0.0);
    TC(1, 2) = std::complex<double>(-st * ct * cp, 0.0);

    TC(2, 0) = TC(0, 2);
    TC(2, 1) = TC(1, 2);
    TC(2, 2) = std::complex<double>(ct * ct, 0.0);
}

void Mult(const MyMatrix &T1, const MyMatrix &TT, MyMatrix &ret) {

    for (int i = 0; i < T1.GetNrows(); ++i) {
        for (int j = 0; j < T1.GetNcols(); ++j) {
            ret(i, j) = std::complex<double>(0.0, 0.0);
            for (int icol = 0; icol < T1.GetNcols(); ++icol) ret(i, j) += (T1(i, icol) * TT(icol, j));
        }
    }
}

void Mult(const MyMatrix &T1, const MyMatrix &TT, const MyMatrix &TTT, MyMatrix &ret) {

    MyMatrix temp;
    Mult(TT, TTT, temp);
    Mult(T1, temp, ret);
}

void Mult(const MyMatrix &T1, const std::complex<double> alpha, MyMatrix &ret) {

    for (int i = 0; i < T1.GetNrows(); ++i) {
        for (int j = 0; j < T1.GetNcols(); ++j) {
            ret(i, j) = T1(i, j) * alpha;
        }
    }
}

void Equal(const MyMatrix &T1, MyMatrix &ret) {

    for (int i = 0; i < T1.GetNrows(); ++i) {
        for (int j = 0; j < T1.GetNcols(); ++j) {
            ret(i, j) = T1(i, j);
        }
    }
}

void Add(const MyMatrix &T1, const MyMatrix &TT, const MyMatrix &TTT, MyMatrix &ret) {

    for (int i = 0; i < T1.GetNrows(); ++i) {
        for (int j = 0; j < T1.GetNcols(); ++j) {
            ret(i, j) = (T1(i, j) + TT(i, j) + TTT(i, j));
        }
    }
}

void Subtract(const MyMatrix &T1, const MyMatrix &TT, MyMatrix &ret) {

    for (int i = 0; i < T1.GetNrows(); ++i) {
        for (int j = 0; j < T1.GetNcols(); ++j) {
            ret(i, j) = (T1(i, j) - TT(i, j));
        }
    }
}

void Conjugate(const MyMatrix &T1, MyMatrix &ret) {

    for (int i = 0; i < T1.GetNrows(); ++i) {
        for (int j = 0; j < T1.GetNcols(); ++j) {
            ret(i, j) = conj(T1(j, i));
        }
    }
}
