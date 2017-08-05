#include "magnfield.h"

std::vector<double> ConstantField::GetB(const double& x_, 
					const double& y_, 
					const double& z_)
{
  std::vector<double> Bret;

  Bret.push_back( BConstant/std::sqrt(3.0) );   
  Bret.push_back( .0/*BConstant/std::sqrt(3.0)*/ ); 
  Bret.push_back( .0/*BConstant/std::sqrt(3.0)*/ ); 
  
  return Bret;
}

std::vector<double> MagneticField::GetBperp(const double& x_,
					    const double& y_, 
					    const double& z_) 
{
  const double x=x_;
  const double y=y_;
  const double z=z_;
  
  const double dist = std::sqrt((x-xSun)*(x-xSun)+(y-ySun)*(y-ySun)+(z-zSun)*(z-zSun));
  
  std::vector<double> Bperp(3,0.0);
  std::vector<double> Btotal = GetB(x,y,z);
  std::vector<double> Versor(3,0.0);
  
  Versor[0] = (x-xSun)/dist;
  Versor[1] = (y-ySun)/dist;
  Versor[2] = (z-zSun)/dist;
  
  Bperp[0] = -(Btotal[1]*Versor[2]-Btotal[2]*Versor[1]);
  Bperp[1] =  (Btotal[0]*Versor[2]-Btotal[2]*Versor[0]);
  Bperp[2] = -(Btotal[0]*Versor[1]-Btotal[1]*Versor[0]);

#ifdef DEBUGMODE  
  const double BtotalDotVersor = Btotal[0]*Versor[0]+Btotal[1]*Versor[1]+Btotal[2]*Versor[2];
  const double BperpDotVersor = Bperp[0]*Versor[0]+Bperp[1]*Versor[1]+Bperp[2]*Versor[2];
  const double BtotalNorm = normVector(Btotal);
  const double BperpNorm = normVector(Bperp);
  const double VersorNorm = normVector(Versor);
  
  std::cout<<std::scientific;
  std::cout<<Versor[0]<<"\t"<<Versor[1]<<"\t"<<Versor[2]<<"\t"<<Bperp[0]<<"\t"<<Bperp[1]<<"\t"<<Bperp[2]<<"\t";
  std::cout<<BtotalDotVersor/BtotalNorm<<"\t"<<std::sin(std::acos(BtotalDotVersor/BtotalNorm))<<"\t"<<BperpNorm/BtotalNorm<<std::endl;
#endif  

  return Bperp;
}


