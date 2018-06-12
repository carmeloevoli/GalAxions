#include "Pshirkov.h"
#include "constants.h"

PshirkovField::PshirkovField(const unsigned int& bMode_, const double& B0Turb_, const double& rscaleTurb_,
        const double& zscaleTurb_)
        :MagneticField() {

    bMode = bMode_;

    p = (bMode == ASS) ? deg2rad(-5.0) : deg2rad(-6.0); // [ASS,BSS]
    z0 = 1.0;  // [kpc]
    d = -0.6; // [kpc]
    B0 = 2.0;  // [muG]
    Rc = 5.0;  // [kpc]
    z0H = 1.3; // [kpc]
    R0H = 8.0; // [kpc]
    B0H[0] = (bMode == ASS) ? 2.0 : 4.0; // [muG] [south]
    B0H[1] = 4.0;  // [mug] [north]
    z1H[0] = 0.25; // [kpc] [inner]
    z1H[1] = 0.40; // [kpc] [outer]
}

double PshirkovField::GetBdisk(const double& rho_, const double& theta_, const double& z_) {
    const double theta = theta_ + M_PI; // CHECK!!

    const double b = 1.0 / tan(p);
    const double phi_disk = b * log(1.0 + d / sun_r) - M_PI / 2.0;

    double Bdisk = cos(theta - b * log(rho_ / sun_r) + phi_disk);

    if (bMode == ASS)
        Bdisk = std::abs(Bdisk);

    Bdisk *= exp(-std::abs(z_) / z0);

    Bdisk *= (rho_ < Rc) ? B0 * sun_r / Rc / cos(phi_disk) : B0 * sun_r / rho_ / cos(phi_disk);

    //Bret[0] = Bdisk * sin(p[mode_]);
    //Bret[1] = Bdisk * cos(p[mode_]) * (-1.);

    return Bdisk;
}

double PshirkovField::GetBhalo(const double& rho_, const double& z_) {
    double z1now = (fabs(z_) < z0H) ? z1H[0] : z1H[1];
    double Bhalo = (z_ < 0.) ? B0H[0] : B0H[1];

    Bhalo /= (1.0 + pow2((fabs(z_) - z0H) / z1now)); // Equation 8
    Bhalo *= rho_ / R0H;
    Bhalo *= exp(1.0 - rho_ / R0H);

    return Bhalo;
}

std::vector<double> PshirkovField::GetB(const double& x_, const double& y_, const double& z_) {
    std::vector<double> Bret(3, 0);

    double Bdisk, Bhalo;

    const double rho = (sqrt(x_ * x_ + y_ * y_) > 1e-3) ? sqrt(x_ * x_ + y_ * y_) : 1e-3;
    const double theta = atan2(y_, x_); // +M_PI/2.0;
    const double z = z_;

    Bdisk = GetBdisk(rho, theta, z);
    Bhalo = GetBhalo(rho, z);

    Bret[0] = (Bdisk + Bhalo) * cos(theta);
    Bret[1] = (Bdisk + Bhalo) * sin(theta);
    Bret[2] = 0.;

    return Bret;
}
