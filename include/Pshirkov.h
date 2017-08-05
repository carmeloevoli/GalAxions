#ifndef _PSHIRKOVFIELD_H
#define _PSHIRKOVFIELD_H

#include <cmath>
#include <vector>
#include <cstdlib>
#include "magnfield.h"
#include "constants.h"

#define ASS 0
#define BSS 1

class PshirkovField : public MagneticField {
  
 private:
  unsigned int bMode; // [ASS,BSS]
  double p; // [ASS,BSS]
  double z0;
  double d;
  double B0;
  double Rc;
  double z0H;
  double R0H;
  double B0H[2]; // [mug] [ASS south, ASS north+BSS]
  double z1H[2]; // [kpc] [inner,outer]
  
 public :
  PshirkovField() { }
  PshirkovField(const unsigned int&,const double&,const double&,const double&);
  virtual ~PshirkovField() { }
  
  virtual std::vector<double> GetB(const double&,const double&,const double&);
  
  double GetBdisk(const double&,const double&,const double&);
  double GetBhalo(const double&,const double&);
};

#endif
