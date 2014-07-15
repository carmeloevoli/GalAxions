#ifndef __FeLine_farrar_h
#define __FeLine_farrar_h

#include <cmath>
#include "magnfield.h"

class FarrarField : public MagneticField {
  
 public :
  FarrarField() { }
  FarrarField(int);
  virtual ~FarrarField() { }
  
  virtual std::vector<double> GetB(double x, double y, double z);
  
 protected :
  double bj[8];
  double fj[8];
  double rj[8];
  double bring;
  double hdisk;
  double wdisk;
  
  double Bn;
  double Bs;
  double rn;
  double rs;
  double wh;
  double z0;
  
  double BX;
  double Theta0X;
  double rcX;
  double rX;
  
  double p;
};

#endif
