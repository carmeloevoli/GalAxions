#ifndef INCLUDE_JF12_H_
#define INCLUDE_JF12_H_

#include <cmath>
#include <limits>
#include "magneticfield.h"

class JF12Field : public MagneticField {
private:
    bool useRegular;
    bool useStriated;
    bool useTurbulent;

    // disk spiral arms
    double rArms[8];       // radii where each arm crosses the negative x-axis
    double pitch;          // pitch angle
    double sinPitch, cosPitch, tan90MinusPitch;

    // Regular field ----------------------------------------------------------
    // disk
    double bDisk[8];       // field strengths of arms at r=5 kpc
    double bRing;          // ring field strength 3<r<5 kpc
    double hDisk, wDisk;   // disk/halo transistion and width
    // toroidal halo
    double bNorth, bSouth; // northern, southern halo field strength
    double rNorth, rSouth; // northern, southern transistion radius
    double wHalo, z0;      // transistion width and vertical scale height
    // poloidal halo
    double bX;             // field strength at origin
    double thetaX0;        // constant elevation angle at r > rXc, z = 0
    double sinThetaX0, cosThetaX0, tanThetaX0;
    double rXc;            // radius of varying elevation angle region
    double rX;             // exponential scale height

    // Striated field ---------------------------------------------------------
    double sqrtbeta;       // relative strength of striated field

    // Turbulent field --------------------------------------------------------
    // disk
    double bDiskTurb[8]; // field strengths in arms at r=5 kpc
    double bDiskTurb5;   // field strength at r<5kpc
    double zDiskTurb;     // Gaussian scale height of disk
    // halo
    double bHaloTurb; // halo field strength
    double rHaloTurb; // exponential scale length
    double zHaloTurb; // Gaussian scale height

public:
    JF12Field();

    void setUseRegular(bool use);

    void setUseStriated(bool use);

    void setUseTurbulent(bool use);

    bool isUsingRegular();

    bool isUsingStriated();

    bool isUsingTurbulent();

    // Regular field component
    std::vector<double> getRegularField(const double& x, const double& y, const double& z) const;

    // Regular and striated field component
    //Vector3d getStriatedField(const Vector3d& pos) const;

    // Brms of the turbulent field
    double getTurbulentStrength(const double& x, const double& y, const double& z) const;

    // Turbulent field component
    std::vector<double> getTurbulentField(const double& x, const double& y, const double& z) const;

    // All set field components
    std::vector<double> GetB(const double& x, const double& y, const double& z) override;
};

#endif /* INCLUDE_JF12_H_ */
