#include "gas.h"

double Cordes91::get(const double& x, const double& y, const double& z) const {
    double r = std::sqrt(x * x + y * y);
    double ne1 = fne1 * exp(-fabs(z) / H1) * exp(-pow2(r / A1));
    double ne2 = fne2 * exp(-fabs(z) / H2) * exp(-pow2((r - R2) / A2));
    return ne1 + ne2;
}

double YMW16::get(const double& x, const double& y, const double& z) const {
    return ymw16_ne(x / pc, y / pc, z / pc) / cm3;
}
