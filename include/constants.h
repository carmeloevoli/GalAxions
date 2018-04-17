#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#include "mks.h"

// ANGLE CONSTANTS

const double PI = M_PI;
const double PiOver2 = PI / 2.0;
const double PiOver4 = PI / 4.0;
const double DegToRad = PI / 180.0;
const double RadToDeg = 180.0 / PI;

// LOS CONSTANTS

const double sun_r = 8.5 * kpc;
const double sun_x = -8.5 * kpc;
const double sun_y = 0;
const double sun_z = 0.;
const double step_size = 0.001 * kpc;
const double x_max = 30 * kpc;
const double y_max = 30 * kpc;
const double z_max = 30. * kpc;

// MAP CONSTANTS

const double maxlong = 180.0;
const double minlong = -180.0;
const double maxlat = 89.75;
const double minlat = -89.75;
const int nlong = 361;
const int nlat = 360;

#endif
