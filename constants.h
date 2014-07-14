#ifndef __CONSTANTS_H
#define __CONSTANTS_H

const double PI =  3.14159265358979323846; 
const double PiOver2 = PI/2.0;
const double PiOver4 = PI/4.0;
const double DegToRad = PI/180.0;
const double RadToDeg = 180.0/PI;
const double kpc2cm = 3.085677580e21;

const bool iHealpix = true;
const int resolution = 7;

// PHYSICAL CONSTANTS ******************************************************************

const double m_e = 5.11; // [100 keV]
const double alpha_em = 1.0/137.0;
const double sigma_Thomson = 6.65e-25;
const double convfactor = 3.06e24;  // cm^-1 --> Mpc^-1

// CONSTANTS FOR THE GALAXY ************************************************************

const double rSun = 8.5;          // Sun radial position
const double xSun = -8.5;
const double ySun = 0;
const double zSun = 0.;           // Sun height
const double ds   = 0.001;        // Step size (kpc)
const double dscm = ds * kpc2cm;  // Step size (cm)
const double xmin = -30;
const double xmax = 30;
const double ymin = -30;
const double ymax = 30;
const double zmin = -30;
const double zmax = 30;

// CONSTANTS FOR THE nuplot GRID ******************************************************

const double maxlong = 180.0;
const double minlong = -180.0;
const double maxlat = 89.75;
const double minlat = -89.75;
const int nlong = 361;
const int nlat = 360;

// *************************************************************************************

inline double sgn(double a) { return (a>=0)-(a<0); }

#endif
