#ifndef _PSHIRKOVFIELD_H
#define _PSHIRKOVFIELD_H

#include <cmath>
#include <vector>
#include "magnfield.h"
#include "constants.h"

class PshirkovField : public MagneticField {
  
 public :
  PshirkovField() { }
  PshirkovField(int, double, double, double, double, double, double, double, double, double, double);
  virtual ~PshirkovField() { }
  
  virtual std::vector<double> GetB(double x, double y, double z);
};

#endif
