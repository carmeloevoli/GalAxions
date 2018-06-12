#ifndef GOLCONDA_MYMATRIX_H_
#define GOLCONDA_MYMATRIX_H_

#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>

class MyMatrix {
protected:
    std::vector<std::vector<std::complex<double> > > elements;

public:
    MyMatrix(bool diag = false) {
        std::vector<std::complex<double> > temp(3, std::complex<double>(0.0, 0.0));
        elements.push_back(temp);
        elements.push_back(temp);
        elements.push_back(temp);

        if (diag) {
            for (int i = 0; i < 3; i++) elements[i][i] = std::complex<double>(1.0, 0.0);
        }
    }

    MyMatrix(const MyMatrix &m) {
        elements = m.elements;
    }

    MyMatrix &operator=(const MyMatrix &m) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                elements[i][j] = m(i, j);
            }
        }
        return *this;
    }

    MyMatrix &operator+=(const MyMatrix &m) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                elements[i][j] += m(i, j);
            }
        }
        return *this;
    }

    ~MyMatrix() { elements.clear(); }

    std::complex<double> operator()(int i, int j) const { return elements[i][j]; }

    std::complex<double> &operator()(int i, int j) { return elements[i][j]; }

    inline int GetNcols() const {
        return 3;
    }

    inline int GetNrows() const {
        return 3;
    }

    void Print() const {
        std::cout << "3x3 Matrix reads : \n" << std::endl;
        std::cout << "------------------ \n" << std::endl;
        std::cout << elements[0][0] << "\t" << elements[0][1] << "\t" << elements[0][2] << std::endl;
        std::cout << elements[1][0] << "\t" << elements[1][1] << "\t" << elements[1][2] << std::endl;
        std::cout << elements[2][0] << "\t" << elements[2][1] << "\t" << elements[2][2] << "\n" << std::endl;
    }
};

void InitD(const double DeltaPerp, const double DeltaPar, const double DeltaA, const double DeltaAg, double D[]);

void InitTA(const double cp, const double sp, MyMatrix &TA);

void InitTB(const double cp, const double sp, const double theta, MyMatrix &TB);

void InitTC(const double cp, const double sp, const double theta, MyMatrix &TC);

void InitTABC(const double cp, const double sp, const double theta, MyMatrix &TA, MyMatrix &TB, MyMatrix &TC);

void Add(const MyMatrix &T, const MyMatrix &TT, const MyMatrix &TTT, MyMatrix &ret);

void Subtract(const MyMatrix &T, const MyMatrix &TT, MyMatrix &ret);

void Mult(const MyMatrix &T, const MyMatrix &TT, const MyMatrix &TTT, MyMatrix &ret);

void Mult(const MyMatrix &T, const MyMatrix &TT, MyMatrix &ret);

void Mult(const MyMatrix &T, const std::complex<double> alpha, MyMatrix &ret);

void Equal(const MyMatrix &T, MyMatrix &ret);

void Conjugate(const MyMatrix &T, MyMatrix &ret);

#endif
