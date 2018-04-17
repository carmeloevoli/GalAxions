#ifndef __FeLine_farrar_h
#define __FeLine_farrar_h

#include <cmath>
#include "magneticfield.h"

class FarrarField : public MagneticField {
  
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

 public :
  FarrarField();// { }
  virtual ~FarrarField() { }
  
  virtual std::vector<double> GetB(const double&, const double&,const double&);
};

#endif
