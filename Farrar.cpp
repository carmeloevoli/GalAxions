#include "Farrar.h"

FarrarField::FarrarField() : MagneticField() {
  
  bj[0] = 0.1;
  bj[1] = 3.0;
  bj[2] = -0.9;
  bj[3] = -0.8;
  bj[4] = -2.0;
  bj[5] = -4.2;
  bj[6] = 0.0;
  bj[7] = 0.0; // 2.7; WHY?
  
  fj[0] = 0.130;
  fj[1] = 0.165;
  fj[2] = 0.094;
  fj[3] = 0.122;
  fj[4] = 0.13;
  fj[5] = 0.118;
  fj[6] = 0.084;
  fj[7] = 0.156;
  
  for (int i = 0; i < 7; i++) bj[7] -= (fj[i]*bj[i]/fj[7]);
  
  rj[0] = 5.1;
  rj[1] = 6.3;
  rj[2] = 7.1;
  rj[3] = 8.3;
  rj[4] = 9.8;
  rj[5] = 11.4;
  rj[6] = 12.7;
  rj[7] = 15.5;
  
  bring = 0.1;
  hdisk = 0.40;
  wdisk = 0.27;
  
  Bn = 1.4;
  Bs = -1.1;
  rn = 9.22;
  rs = 16.7;
  wh = 0.2;
  z0 = 5.3;
  
  BX = 4.6;
  Theta0X = 49.0*DegToRad;
  rcX = 4.8;
  rX = 2.9;
  
  p = 11.5*DegToRad;
}

std::vector<double> FarrarField::GetB(const double& x_, 
				      const double& y_, 
				      const double& z_)
{
  std::vector<double> Bret(3,0);
  
  const double r = std::sqrt(x_*x_+y_*y_);
  if (r > 20 || std::sqrt(r*r+z_*z_) < .1) return Bret;
  
  double Bdisk=0;
  double Bhalo=0;
  double Bxx=0;
  
  double rp=0;
  double ThetaX=0;
  double phi = std::atan2(y_,x_);// + M_PI/2.0;
  double sinphi = sin(phi);
  double cosphi = cos(phi);
  double Ldisk = 1.0/(1.0+exp(-2.0*(fabs(z_)-hdisk)/wdisk));
  
  if ( r > 3.0 ){
    if ( r < 5.0 ){
      Bret[0] = -bring*sinphi*(1.0-Ldisk);
      Bret[1] = bring*cosphi*(1.0-Ldisk);
    }
    else {
      const double tanfactor = 1.0/tan(PiOver2-p);
      double rxneg = r*exp(-(phi-PI)*tanfactor);
      if (rxneg > rj[7]) rxneg = r*exp(-(phi+PI)*tanfactor);   
      if (rxneg > rj[7]) rxneg = r*exp(-(phi+3.0*PI)*tanfactor);   

      for (int loopcounter = 7; loopcounter>=0; loopcounter--) { 
	if (rxneg < rj[loopcounter]) Bdisk = bj[loopcounter]; }
      
      Bdisk *= (5.0/r) ;
      Bret[0] = Bdisk*sin(p-phi)*( 1.0 - Ldisk );
      Bret[1] = Bdisk*cos(p-phi)*( 1.0 - Ldisk );
    }
  }
  
  double Lhalo = 0;
  
  if (z_>=0) { 
    Bhalo = Bn;
    Lhalo = 1.0 / (1.0 + exp(-2.0*(r-rn)/wh) );
  }
  else { 
    Bhalo = Bs; 
    Lhalo = 1.0 / (1.0 + exp(-2.0*(r-rs)/wh) ); 
  }
  
  Bhalo *= (exp(-fabs(z_)/z0)*Ldisk*(1.0-Lhalo));
  Bret[0] -= (Bhalo*sinphi);
  Bret[1] += (Bhalo*cosphi);
  
  double ztheta0x = fabs(z_)/tan(Theta0X);
  double rpcX = rcX+ztheta0x;
  
  if (r<rpcX){ // interior region, with varying elevation angle
    rp   = r*rcX/rpcX;
    ThetaX = (z_!=0) ? atan(fabs(z_)/(r-rp)) : PiOver2;
    Bxx = BX*exp(-rp/rX)*pow(rcX/rpcX,2);
  }
  else { // exterior region with constant elevation angle
    rp = r-ztheta0x;
    ThetaX = Theta0X;
    Bxx = BX*exp(-rp/rX)*(rp/r);
  }
  
  if (z_<0) {
    Bret[0] -= (Bxx*cosphi*cos(ThetaX));
    Bret[1] -= (Bxx*sinphi*cos(ThetaX));
  }
  else {
    Bret[0] += (Bxx*cosphi*cos(ThetaX));
    Bret[1] += (Bxx*sinphi*cos(ThetaX));
  }
  Bret[2] += (Bxx*sin(ThetaX));
  
  //cout<<fixed<<setprecision(2)<<x<<"\t"<<y<<"\t"<<z<<"\t"<<Bret[0]<<"\t"<<Bret[1]<<"\t"<<Bret[2]<<endl;
  return Bret;
}


