#ifndef __MAGNFIELD_H
#define __MAGNFIELD_H

#include <iostream>
#include <map>
#include <vector>
#include <cmath>

#include "constants.h"

class MagneticField {
  
 protected :
  int Dim;
  std::vector<double> x_grid;
  std::vector<double> y_grid;
  std::vector<double> z_grid;
  std::map<int, std::vector< std::vector< std::vector<double> > > > Breg;
  std::vector< std::vector< std::vector<double> > > Brand;
      
 public :
  MagneticField() { }
  MagneticField(int N) { 
    Dim = N;
    x_grid = std::vector<double>(N,0.0);
    y_grid = std::vector<double>(N,0.0);
    z_grid = std::vector<double>(N,0.0);
    
    for (int i = 0; i < N; i++) {
      x_grid[i] = xmin + (xmax-xmin)*double(i)/double(N-1);
      y_grid[i] = ymin + (ymax-ymin)*double(i)/double(N-1);
      z_grid[i] = zmin + (zmax-zmin)*double(i)/double(N-1);
    }
    
    for (int i = 0; i < 3; i++) 
      Breg[i]  = std::vector< std::vector< std::vector<double> > >(N, std::vector< std::vector<double> >(N, std::vector<double>(N,0.0)));
    Brand = std::vector< std::vector< std::vector<double> > >(N, std::vector< std::vector<double> >(N, std::vector<double>(N,0.0)));
  }
  
  virtual ~MagneticField() { }
  
  int GetDim() const { return Dim; }
  std::map<int, std::vector< std::vector< std::vector<double> > > >& GetB() { return Breg; }
  std::vector<double> GetBperp(double x, double z, double y);
  virtual std::vector<double> GetB(double x, double y, double z) { return std::vector<double>(3,0); }
  double GetBrand(double x, double z, double y);
  void Add(MagneticField * B);
};

class ConstantField : public MagneticField {
  
 public:
  ConstantField() { }
 ConstantField(int N) :
  MagneticField(N) {
    for (unsigned int i = 0; i < N; i++) {
      for (unsigned int j = 0; j < N; j++) {
	for (unsigned int k = 0; k < N; k++) {
	  Breg[0][i][j][k] = 1.0; // muG
	  Breg[1][i][j][k] = 1.0; // muG
	  Breg[2][i][j][k] = 1.0; // muG
	}
      }
    }
  }
  
  ~ConstantField() { }
};

#endif
