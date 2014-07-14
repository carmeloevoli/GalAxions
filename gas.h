#ifndef __GAS_DISTRIBUTION_H
#define __GAS_DISTRIBUTION_H

#include <cmath>
#include "constants.h"

class Ferriere {
  
 private:
  int model;
  
 public:
  Ferriere() { model=0; }
  ~Ferriere() { }
  
  int getModel(){ 
    return model; 
  };
  
  double getElectrons(const double& rGal, const double& zGal){ 
    return 2.13e-2/1e-7*std::exp( -(rGal-rSun)/18.24 ) * std::min( 1.0, std::exp( (fabs(zGal)-1.0)/0.5) ); 
  } 
  
  double getH2(const double& rGal, const double& zGal){ 
    return 2.0*(4.06*exp(-rGal/2.57-fabs(zGal)/0.08)) ;
  }
  
  double getHI(const double& rGal, const double& zGal){ 
    if (rGal >= 2.75) 
      return 0.32*exp(-rGal/18.24 - fabs(zGal)/0.52);
    else 
      return 0.0;
  }
};

#endif
