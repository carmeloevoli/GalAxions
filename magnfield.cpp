#include "magnfield.h"

void MagneticField::Add(MagneticField* B) {
  
  if (Dim != B->GetDim()) {
    std::cerr<<"Trying to sum incompatible magnetic fields!"<<std::endl;
    return ;
  }
  
  std::map<int, std::vector< std::vector< std::vector<double> > > > Br1 = B->GetB();
  
  for (int dir = 0; dir < 3; dir++) {
    for (int i = 0; i < Dim; i++) {
      for (int j = 0; j < Dim; j++) {
	for (int k = 0; k < Dim; k++) {
	  Breg[dir][i][j][k] += Br1[dir][i][j][k];
	}
      }
    }
  }
  
  return ;
}

std::vector<double> MagneticField::GetBperp(double x, double z, double y) {
  
  std::vector<double> Bperp(3,0.0);
  std::vector<double> Binterp = GetB(x,y,z);
  
  double dist = std::sqrt((x-xSun)*(x-xSun)+(y-ySun)*(y-ySun)+(z-zSun)*(z-zSun));

  Bperp[0] = -(Binterp[1]*(z)-Binterp[2]*(y))/dist*1e3;
  Bperp[1] = (Binterp[0]*(z)-Binterp[2]*(x-xSun))/dist*1e3;
  Bperp[2] = -(Binterp[0]*(y)-Binterp[1]*(x-xSun))/dist*1e3;

  return Bperp;
}

double MagneticField::GetBrand(double x, double z, double y) {
  
  unsigned int i = int(floor((x-x_grid[0])/(x_grid[1]-x_grid[0])));
  unsigned int j = int(floor((y-y_grid[0])/(y_grid[1]-y_grid[0])));
  unsigned int k = int(floor((z-z_grid[0])/(z_grid[1]-z_grid[0])));
  if (i > x_grid.size()-2) i = x_grid.size()-2;
  if (j > y_grid.size()-2) j = y_grid.size()-2;
  if (k > z_grid.size()-2) k = z_grid.size()-2;
  
  const double t = (x-x_grid[i])/(x_grid[i+1]-x_grid[i]);
  const double u = (y-y_grid[j])/(y_grid[j+1]-y_grid[j]);
  const double v = (z-z_grid[k])/(z_grid[k+1]-z_grid[k]);
  
  double Bfield_low = Brand[i][j][k]*(1.0-t)*(1.0-u) + Brand[i+1][j][k]*(t)*(1.0-u) + Brand[i][j+1][k]*(1.0-t)*(u) + Brand[i+1][j+1][k]*(t)*(u);
  
  double Bfield_high = Brand[i][j][k+1]*(1.0-t)*(1.0-u) + Brand[i+1][j][k+1]*(t)*(1.0-u) + Brand[i][j+1][k+1]*(1.0-t)*(u) + Brand[i+1][j+1][k+1]*(t)*(u);
  
  return Bfield_low*(1.0-v)+Bfield_high*v;
}
