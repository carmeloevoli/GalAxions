#ifndef INCLUDE_CGS_H_
#define INCLUDE_CGS_H_

#include <cmath>

#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define pow4(A) ((A)*(A)*(A)*(A))
#define pow5(A) ((A)*(A)*(A)*(A)*(A))
#define sech(A) (1.0 / std::cosh(A))
#define sech2(A) (1.0 / std::cosh(A) / std::cosh(A))
#define deg2rad(A) (A * M_PI / 180.0)

// CGS units
static const double centimeter = 1;
static const double gram = 1;
static const double second = 1;
static const double gauss = 1;
static const double erg = 1;

// length units
static const double meter = 1e2 * centimeter;
static const double parsec = 3.0856775807e18 * centimeter;
static const double kiloparsec = 1e3 * parsec;
static const double megaparsec = 1e6 * parsec;
static const double gigaparsec = 1e9 * parsec;

// magnetic field units
static const double microgauss = 1e-6 * gauss;
static const double nanogauss = 1e-9 * gauss;

// energy
static const double joule = 1e7 * erg;
static const double electronvolt = 1.60217657e-19 * joule;
static const double kiloelectronvolt = 1e3 * electronvolt;
static const double megaelectronvolt = 1e6 * electronvolt;
static const double gigaelectronvolt = 1e9 * electronvolt;
static const double teraelectronvolt = 1e12 * electronvolt;
static const double petaelectronvolt = 1e15 * electronvolt;
static const double exaelectronvolt = 1e18 * electronvolt;

// angle
static const double PiOver2 = M_PI / 2.0;
static const double PiOver4 = M_PI / 4.0;

// abbreviations
static const double pc = parsec;
static const double kpc = kiloparsec;
static const double Mpc = megaparsec;
static const double Gpc = gigaparsec;
static const double cm = centimeter;
static const double cm2 = cm * cm;
static const double cm3 = cm * cm * cm;
static const double muG = microgauss;
static const double nG = nanogauss;
static const double eV = electronvolt;
static const double keV = kiloelectronvolt;
static const double MeV = megaelectronvolt;
static const double GeV = gigaelectronvolt;
static const double TeV = teraelectronvolt;
static const double PeV = petaelectronvolt;

#endif /* INCLUDE_CGS_H_ */
