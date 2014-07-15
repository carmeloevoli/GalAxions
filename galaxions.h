#ifndef __GALAXIONS_H
#define __GALAXIONS_H

#include <iostream>
#include <iomanip>
#include <string>

#include "gas.h"
#include "magnfield.h"
#include "Farrar.h"
#include "Pshirkov.h"

enum MagneticFieldType { CONSTANT, FARRAR, PSHIRKOV };

struct los_struct {
  double l;
  double b;
  int nSteps;
  std::vector<double> magneticField;
  std::vector<double> distance;
  std::vector<double> psik;
  std::vector<double> electronDensity;
  std::vector<double> H2Density;
  std::vector<double> HIDensity;
};

class galAxions{
 private:
  
  double mass; 
  double gag;

  std::string initFilename;
  //std::string rootFilename;
  std::string asciiFilename;
  std::string probabilyMapFilename; // FITS file
  std::string pdMapFilename; // FITS file
  std::string paMapFilename; // FITS file
  
  MagneticField * magneticField; // = NULL;
  Ferriere * gas;

  los_struct los;
  
 public:
  
  galAxions(const double& _mass, 
	    const double& _gag,
	    const std::string& _initFilename)
    { 
      mass = _mass; 
      gag = _gag;
      initFilename = _initFilename;
      asciiFilename = _initFilename; 
      asciiFilename += ".dat";
      probabilyMapFilename = _initFilename;
      probabilyMapFilename += ".fits";
      pdMapFilename = _initFilename;
      pdMapFilename += "_PD.fits"; 
      paMapFilename = _initFilename;
      paMapFilename += "_PA.fits";
      //magnfield = NULL;
    }
  
  ~galAxions()
    {
      if (magneticField) delete magneticField;
      if (gas) delete gas;
    }
  
  double getMass(){ return mass; }
  double getGag(){ return gag; }
  std::string getInitFilename(){ return initFilename; }
  
  void createGasDensity();  
  void createMagneticField(const MagneticFieldType&);
  void createLos(const double&, const double&);
};

#endif
