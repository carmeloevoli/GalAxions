#ifndef __MAGNFIELD_H
#define __MAGNFIELD_H

#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>

#include "constants.h"

#define POW2(A) ((A)*(A))

inline double normVector( const std::vector<double>& v ){ return std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }

class MagneticField {
  
 private :
  double Bmin;
  double Bmax;
  
 public :
  MagneticField() { }
  
  virtual ~MagneticField() { }
  
  std::vector<double> GetBperp(const double&,const double&,const double&);
  
  double GetBtransverse(const double&,const double&,const double&);
  
  virtual std::vector<double> GetB(const double&,const double&,const double&) { return std::vector<double>(3,0); }
  
};

class ConstantField : public MagneticField {
  
 private:
  double BConstant;
  
 public:
  ConstantField() { }
  ConstantField(const double& B_){ BConstant = B_; }
  virtual ~ConstantField() { }
  
  virtual std::vector<double> GetB(const double&,const double&,const double&);
};

#endif
