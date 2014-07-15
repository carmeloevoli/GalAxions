#include "Pshirkov.h"
#include "constants.h"

PshirkovField::PshirkovField(int N, 
			     double B0, 
			     double RC, 
			     double Z0, 
			     double R0, 
			     double B0H, 
			     double R0H, 
			     double Z0H, 
			     double B0turb, 
			     double rscale_turb, 
			     double zscale_turb) : MagneticField(N) {
  
  double Bdisk, Bhalo;
  
  //double B0 = 2.0, RC = 5.0, Z0 = 1.0, R0 = 10.0; 
  //double B0H = 4., R0H = 8., Z0H = 1.3;
  
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double r = std::sqrt(x_grid[i]*x_grid[i]+y_grid[j]*y_grid[j]);
      if (fabs(r) < 1e-3) r = 1e-3;
      
      for (int k = 0; k < N; k++) {

	double phi = atan2(y_grid[j],x_grid[i])+M_PI/2.0;
	double z = z_grid[k];
        
	double Z1H = (fabs(z) < Z0H) ? 0.2 : 0.4; 
	
	Bdisk = B0 * ((r < RC) ? exp(-fabs(z)/Z0) : exp(-(r-rSun)/R0-fabs(z)/Z0));
	Bhalo = B0H / (1.+pow((fabs(z)-Z0H)/Z1H,2)) * r/R0H * exp(-1.-r/R0H);
	
	Breg[0][i][j][k] = (Bdisk+Bhalo)*cos(phi);
	Breg[1][i][j][k] = (Bdisk+Bhalo)*sin(phi);
	Breg[2][i][j][k] = 0.;
        
	Brand[i][j][k] = B0turb*exp(-(r-rSun)/rscale_turb)*exp(-fabs(z/zscale_turb));
      }
    }
  }
}

std::vector<double> PshirkovField::GetB(double x, double y, double z) {
  
  std::vector<double> Bret(3,0);
  
  double Bdisk, Bhalo;
  
  // Disk parameters from Table 3
  
  const double B0 = 2.0;
  const double RC = 5.0;
  const double Z0 = 1.0;
  const double R0 = 10.0; // WHERE??
  
  // Halo parameters from Table 3
  
  const double B0H = 4.0;
  const double Z0H = 1.3;
  const double Z1H = (fabs(z) < Z0H) ? 0.25 : 0.4; 
  const double R0H = 8.0;
  
  double r = std::sqrt(x*x+y*y);
  if (fabs(r) < 1e-3) r = 1e-3;
  
  double phi = atan2(y,x)+M_PI/2.0;
  
  Bdisk = B0 * ( (r <= RC) ? exp(-fabs(z)/Z0) : exp(-(r-rSun)/R0-fabs(z)/Z0) ); // Eq. 6 
  Bhalo = B0H / (1.0+pow((fabs(z)-Z0H)/Z1H, 2)) * r/R0H * exp(1.0-r/R0H); // Eq. 8
  
  Bret[0] = (Bdisk+Bhalo)*cos(phi);
  Bret[1] = (Bdisk+Bhalo)*sin(phi);
  Bret[2] = 0.;
  
  return Bret;  
}
