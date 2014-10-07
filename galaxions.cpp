#include "galaxions.h"

#define BASSAN

void galAxions::createMagneticField(const MagneticFieldType& btype_,
				    const int& bmode_)
{
  if ( btype_ == CONSTANT ){
    std::cout<<"#Init constant magnetic field!"<<std::endl;
    magneticField = new ConstantField(2.0);
  }
  else if ( btype_ == FARRAR ){
    std::cout<<"#Init Farrar magnetic field!"<<std::endl;
    magneticField = new FarrarField();
  }
  else if ( btype_ == PSHIRKOV ){
    std::cout<<"#Init Pshirkov magnetic field!"<<std::endl;
    magneticField = new PshirkovField(bmode_, 0.0, 8.5, 2.0);
  }

  return;
}

void galAxions::createGasDensity(void)
{
  std::cout<<"#Init gas!"<<std::endl;
  gas = new Ferriere();
  
  return;
}

void galAxions::printLos(const double& rmax_)
{
  galaxyStream<<std::scientific<<std::setprecision(5);
  
  for ( unsigned int i=0;i<los.nSteps;i++ ){
    if ( los.distance[i] < rmax_ ){
      galaxyStream<<los.distance[i]<<"\t";
      galaxyStream<<los.magneticFieldPerp[i]<<"\t";
      galaxyStream<<los.magneticFieldTotal[i]<<"\t";
      galaxyStream<<los.psik[i]<<"\t";
      galaxyStream<<los.electronDensity[i]<<"\t";
      galaxyStream<<los.H2Density[i]<<"\t";
      galaxyStream<<los.HIDensity[i]<<std::endl;
    }
  }
  
  return;
}

void galAxions::createLos(const double& ldeg_, 
			  const double& bdeg_)
{ 
  los.ldeg = ldeg_;
  los.bdeg = bdeg_;
  
  los.l = los.ldeg*DegToRad; // \phi in [0,2\pi] 
  //los.b = (bdeg < 0) ? (90.0-bdeg)*DegToRad : -1; // \theta in [0,\pi]
  los.b = los.bdeg*DegToRad;
    
  double sinl = sin( los.l );
  double cosl = cos( los.l );
   
  double sinb = sin( los.b );
  double cosb = cos( los.b );
  
  double cosbcosl = cosb*cosl;
  double cosbsinl = cosb*sinl;
    
  double distanceAlongLos = 0.0;
  double xGalactoCentric = 0.0;
  double yGalactoCentric = 0.0;
  double zGalactoCentric = 0.0;
  
  std::vector<double> ydir;
  /*std::vector<double> xdir(3,0.0);
  std::vector<double> dir(3,0.0);
  dir[0] = -cosbcosl;
  dir[1] = -cosbsinl;
  dir[2] = -sinb;*/
  
  bool done = false;
  
  std::vector<double> Bperp;
  std::vector<double> Btotal, Btest;
  
  while ( fabs(xGalactoCentric) < xmax && fabs(yGalactoCentric) < ymax && fabs(zGalactoCentric) < zmax ){
    
    distanceAlongLos += ds; // [kpc]
    
    if (xSun<0) 
      xGalactoCentric = xSun + distanceAlongLos*cosbcosl;
    else  
      xGalactoCentric = xSun - distanceAlongLos*cosbcosl;
    
    yGalactoCentric = distanceAlongLos*cosbsinl;
    
    zGalactoCentric = distanceAlongLos*sinb;
    
    Bperp.clear();
    Bperp = magneticField->GetBperp(xGalactoCentric,yGalactoCentric,zGalactoCentric); // [muG]
    Btotal = magneticField->GetB(xGalactoCentric,yGalactoCentric,zGalactoCentric); // [muG]

#ifdef FARRARPLOT    
    for ( double ztemp=-3;ztemp<=3;ztemp+=0.01 ){
      Btest = magneticField->GetB(-8.5,0.0,ztemp);
      std::cout<<ztemp<<"\t"<<normVector(Btest)<<"\t";
      Btest = magneticField->GetB(-10.0,0.0,ztemp);
      std::cout<<normVector(Btest)<<"\t"<<std::endl;
    }
    exit(1);
#endif

    los.magneticFieldPerp.push_back( normVector(Bperp) );
    los.magneticFieldTotal.push_back( normVector(Btotal) );
    
    double cosarg = 0;
    double testpsik = 0; // the angle between Btransverse and y axis
    
    if ( !done ){ // Fix once the reference direction
      
      done = true;
      
      ydir = Bperp;
      ydir[0] /= los.magneticFieldPerp.back();
      ydir[1] /= los.magneticFieldPerp.back();
      ydir[2] /= los.magneticFieldPerp.back();
      
      /*xdir[0] = (ydir[1]*dir[2]-dir[1]*ydir[2]);
      xdir[1] = -(ydir[0]*dir[2]-dir[0]*ydir[2]);
      xdir[2] = (ydir[0]*dir[1]-dir[0]*ydir[1]);*/
    }
    
    if ( los.magneticFieldPerp.back() != 0 ){
      cosarg = (ydir[0]*Bperp[0]+ydir[1]*Bperp[1]+ydir[2]*Bperp[2])/los.magneticFieldPerp.back();
      if (cosarg >= 0.99999999999) testpsik = 0.0;
      else if (cosarg <= -0.99999999999) testpsik = PI;
      else testpsik = acos(cosarg);
    }
    else testpsik = 0.0;
    
    los.psik.push_back(testpsik);
    los.distance.push_back(distanceAlongLos);  // [kpc]
    
    double rGalactoCentric = std::sqrt(xGalactoCentric*xGalactoCentric+yGalactoCentric*yGalactoCentric);
    
    los.electronDensity.push_back( gas->getElectrons(rGalactoCentric,zGalactoCentric) );
    los.H2Density.push_back( gas->getH2(rGalactoCentric,zGalactoCentric) );
    los.HIDensity.push_back( gas->getHI(rGalactoCentric,zGalactoCentric) );
  }
  
  los.nSteps = los.distance.size();

  printLos(10.0);
 
  std::reverse(los.distance.begin(),los.distance.end());
  std::reverse(los.magneticFieldPerp.begin(),los.magneticFieldPerp.end());
  std::reverse(los.magneticFieldTotal.begin(),los.magneticFieldTotal.end());
  std::reverse(los.psik.begin(),los.psik.end());
  std::reverse(los.electronDensity.begin(),los.electronDensity.end());
  std::reverse(los.H2Density.begin(),los.H2Density.end());
  std::reverse(los.HIDensity.begin(),los.HIDensity.end());
  
  return;
}

void galAxions::calculateProbability(const unsigned int& nEnergy_,
				     const double& EminInEv_, // [eV]
				     const double& EmaxInEv_, // [eV]
				     const bool& WithDamping_) // [eV]
{
  bool PRINT_SOLVER = true;
  bool PRINT_AT_EACH_TIMESTEP = false; //true;
  
  time_t timeBegin, timeEnd;

  const int nDomains = los.distance.size();
  
  time(&timeBegin);

  for ( unsigned int iE = 0; iE < nEnergy_; iE++ ){
    
    const double EnergyInEv = ( nEnergy_ == 1 ) ? EminInEv_ : pow(10, log10(EminInEv_)+double(iE)/double(nEnergy_-1)*log10(EmaxInEv_/EminInEv_)); 
    
    const double x = EnergyInEv/15.4;
    
    const double XsecH2 = (WithDamping_) ? 45.57e-24*(1.0-2.003*pow(x,-0.5) -4/806/x+50.577*pow(x,-1.5)-171.044*pow(x,-2)+231.608*pow(x,-2.5)-81.885*pow(x,-3))/pow(EnergyInEv/1e3,3.5) : 0.0; // CHECK DIMENSION!
    
    const double XsecHI = (WithDamping_) ? XsecH2/2.833 : 0.0; // CHECK DIMENSION!
    
    const double mass   = axionMassInEv/1e-13;
    const double gag    = gagInGeV/1e-11;
    const double Energy = EnergyInEv/1e5;
    
    const double DeltaA = -7.8e-3*mass*mass/Energy*MpcInverse2kpcInverse; // [kpc^-1] Eq. 3.13
    
    for ( unsigned int i = 0; i < nDomains; i++ ){
      const double Bnow = los.magneticFieldPerp[i]*muG2nG;
      const double nenow = los.electronDensity[i]/1e-7;
      DeltaAgamma.push_back( 1.52e-2*gag*Bnow*MpcInverse2kpcInverse ); // [kpc^-1] Eq. 3.12
      DeltaPl.push_back( -1.1e-4/Energy*nenow*MpcInverse2kpcInverse ); // [kpc^-1] Eq. 3.14
      DeltaQED.push_back( 4.1e-16*Energy*pow(Bnow,2.0)*MpcInverse2kpcInverse ); // [kpc^-1] Eq. 3.15
      DeltaPar.push_back( DeltaPl.back() + 3.5*DeltaQED.back() ); // [kpc^-1] Eq. 3.10
      DeltaPerp.push_back( DeltaPl.back() + 2.0*DeltaQED.back() ); // [kpc^-1] Eq. 3.11
    } // Initialization
    
    MyMatrix OldRho;
    OldRho(0,0) = std::complex<double>(0.0,0.0);
    OldRho(1,1) = std::complex<double>(0.0,0.0); // Photon unpolarized. Initial condition.
    OldRho(2,2) = std::complex<double>(1.0,0.0); // Full ALPs beam. Initial condition.
    
    MyMatrix Tktotal, 
      Rho;
    
#ifdef BASSAN    
    if (PRINT_SOLVER){
      std::cout<<"#Using BassanSolver!"<<std::endl;
      PRINT_SOLVER = false;
    }
    
    BassanSolver(nDomains, 
		 DeltaAgamma, 
		 DeltaPl, 
		 DeltaQED, 
		 DeltaPar, 
		 DeltaPerp, 
		 DeltaA, 
		 los.psik, 
		 los.H2Density, 
		 los.HIDensity, 
		 XsecH2, 
		 XsecHI, 
		 OldRho, 
		 Rho, 
		 Tktotal);
#endif
    
    if ( PRINT_AT_EACH_TIMESTEP ){
      std::cout<<"Energy = "<<EnergyInEv<<std::endl;
      std::cout<<"Tk final "<<std::endl;
      Tktotal.Print();
    }
    
    PagPDF = real(Rho(2,2));
    IavPDF = real(Rho(0,0) + Rho(1,1));
    QavPDF = real(Rho(0,0) - Rho(1,1));
    UavPDF = real(Rho(0,1) + Rho(1,0));
    VavPDF = imag(Rho(1,0) - Rho(0,1));
    
    PolDegPDF = sqrt(pow(QavPDF, 2) + pow(UavPDF, 2) /*+ pow(VavPDF, 2)*/)/IavPDF; // Eq. 3.44
    PosAnglePDF = RadToDeg*0.5*atan2(UavPDF, QavPDF);
    
    if ( IavPDF < 0 ){
      std::cout<<"Warning: Negative Intensity!"<<std::endl;
      Rho.Print();
      std::cout<<"Coordinates: "<<los.ldeg<<" "<<los.bdeg<<std::endl;
    }
    
    outputStream<<std::scientific<<EnergyInEv<<"\t"<<IavPDF<<"\t"<<PagPDF<<std::endl;
    
    DeltaAgamma.clear();
    DeltaPl.clear();
    DeltaQED.clear();
    DeltaPar.clear();
    DeltaPerp.clear();
    
  } //close energy loop

  time(&timeEnd);
  
  std::cout<<"Ended in "<<difftime(timeEnd,timeBegin)<<" seconds."<<std::endl;// in "<<tf-te<<" seconds."<<endl;
  return;
}
