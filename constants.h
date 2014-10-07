#ifndef __CONSTANTS_H
#define __CONSTANTS_H

const double PI =  3.14159265358979323846; 
const double PiOver2 = PI/2.0;
const double PiOver4 = PI/4.0;
const double DegToRad = PI/180.0;
const double RadToDeg = 180.0/PI;

const bool iHealpix = true;
const int resolution = 7;

// PHYSICAL CONVERSIONS ******************************************************************

const double cm2kpc = 1;
const double cmInverse2kpcInverse = 3.06e21;
const double MpcInverse2kpcInverse = 1e-3;
const double kpc2cm = 3.085677580e21;
const double muG2nG = 1e3;

// CONSTANTS FOR THE GALAXY ************************************************************

const double rSun = 8.5;          // [kpc] Sun radial position
const double xSun = -8.5;         // [kpc]
const double ySun = 0;            // [kpc]
const double zSun = 0.;           // [kpc] 
const double ds   = 0.001;        // [kpc] Step size 
const double xmin = -30;          // [kpc] 
const double xmax = 30;           // [kpc] 
const double ymin = -30;          // [kpc] 
const double ymax = 30;           // [kpc] 
const double zmin = -30;          // [kpc] 
const double zmax = 30;           // [kpc] 

// CONSTANTS FOR THE GRID ******************************************************

const double maxlong = 180.0;
const double minlong = -180.0;
const double maxlat = 89.75;
const double minlat = -89.75;
const int nlong = 361;
const int nlat = 360;

// *************************************************************************************

inline double sgn(double a) { return (a>=0)-(a<0); }

#endif
