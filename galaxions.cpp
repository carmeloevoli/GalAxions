#include "galaxions.h"

void galAxions::createMagneticField(const MagneticFieldType& _btype)
{
  if ( _btype == CONSTANT ){
    std::cout<<"#Init constant magnetic field!"<<std::endl;
    magneticField = new ConstantField(100);
  }
  else if ( _btype == FARRAR ){
    std::cout<<"#Init Farrar magnetic field!"<<std::endl;
    magneticField = new FarrarField(0);
  }
  else if ( _btype == PSHIRKOV ){
    std::cout<<"#Init Pshirkov magnetic field!"<<std::endl;
  }

  return;
}

void galAxions::createGasDensity()
{
  std::cout<<"#Init gas!"<<std::endl;
  gas = new Ferriere();
  
  return;
}

void galAxions::createLos(const double& ldeg, const double& bdeg)
{ 
  los.l = ldeg*DegToRad; // \phi in [0,2\pi] 
  los.b = (bdeg < 0) ? (90.0-bdeg)*DegToRad : -1; // \theta in [0,\pi]
    
  //getArrays(los, magneticField, distance_array, Bfield, psik, el_dens, gas_densH2, gas_densHI);
  
  double sinl = sin( los.l );
  double cosl = cos( los.l );
   
  double sinb = sin( los.b );
  double cosb = cos( los.b );
  
  double cosbcosl = cosb*cosl;
  double cosbsinl = cosb*sinl;
    
  double d = 0.0;
  float zz = 0.0;
  float xs = 0.0;
  float ys = 0.0;
  
  std::vector<double> ydir;
  std::vector<double> xdir(3,0.0);
  std::vector<double> dir(3,0.0);
  
  dir[0] = -cosbcosl;
  dir[1] = -cosbsinl;
  dir[2] = -sinb;
  
  bool done = false;
  
  std::vector<double> Bperp;
  
  while (fabs(xs) < xmax && fabs(ys) < ymax && fabs(zz) < zmax) {
    
    d += ds;
    //cout<<d<<"\t"<<ds<<"\t"<<xs<<"\t"<<ys<<"\t"<<zz<<endl;
    
    zz=d*sinb;
    
    if (xSun<0) 
      xs = xSun + d*cosbcosl;
    else  xs = xSun - d*cosbcosl;
    
    ys = d*cosbsinl;
    
    Bperp = magneticField->GetBperp(xs,zz,ys); // [nG]
    
    los.magneticField.push_back( std::sqrt(Bperp[0]*Bperp[0] + Bperp[1]*Bperp[1] + Bperp[2]*Bperp[2]) );
    
    //cerr<<scientific<<d<<"\t"<<Bfield.back()<<endl;
    
    double cosarg = 0;
    double testpsik = 0;
    
    if (!done) { // Fix once the reference direction

      done = true;
      
      ydir = Bperp;
      ydir[0] /= los.magneticField.back();
      ydir[1] /= los.magneticField.back();
      ydir[2] /= los.magneticField.back();
      
      xdir[0] = (ydir[1]*dir[2]-dir[1]*ydir[2]);
      xdir[1] = -(ydir[0]*dir[2]-dir[0]*ydir[2]);
      xdir[2] = (ydir[0]*dir[1]-dir[0]*ydir[1]);
    }
    
    if ( los.magneticField.back() != 0 ){
      cosarg = (ydir[0]*Bperp[0]+ydir[1]*Bperp[1]+ydir[2]*Bperp[2])/los.magneticField.back();
      if (cosarg >= 0.99999999999) testpsik = 0.0;
      else if (cosarg <= -0.99999999999) testpsik = PI;
      else testpsik = acos(cosarg);
    }
    else testpsik = 0.0;
      
    los.psik.push_back(testpsik);
    los.distance.push_back(d*1e-3);  // [Mpc]
    
    //std::cerr<<std::scientific<<los.distance.back()<<"\t"<<los.psik.back()<<std::endl;
    
    double rs = std::sqrt(xs*xs+ys*ys);
    
    los.electronDensity.push_back( gas->getElectrons(rs,zz) );
    los.H2Density.push_back( gas->getH2(rs,zz) );
    los.HIDensity.push_back( gas->getHI(rs,zz) );
  }

  los.nSteps = los.distance.size();

  std::cout<<std::scientific<<std::setprecision(3);

  for ( unsigned int i=0;i<los.nSteps;i++ ){
    std::cout<<los.distance[i]<<"\t";
    std::cout<<los.magneticField[i]<<"\t";
    std::cout<<los.psik[i]<<"\t";
    std::cout<<los.electronDensity[i]<<"\t";
    std::cout<<los.H2Density[i]<<"\t";
    std::cout<<los.HIDensity[i]<<std::endl;
  }

  return;
}
