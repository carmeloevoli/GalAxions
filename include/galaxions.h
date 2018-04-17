#ifndef __GALAXIONS_H
#define __GALAXIONS_H

#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <complex>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gas.h"
#include "magnfield.h"
#include "Farrar.h"
#include "Pshirkov.h"
#include "mymatrix.h"
#include "solvers.h"

enum MagneticFieldType { CONSTANT, FARRAR, PSHIRKOV };

struct los_struct {
  double l;
  double ldeg;
  double b;
  double bdeg;
  int nSteps;
  std::vector<double> magneticFieldPerp;
  std::vector<double> magneticFieldTotal;
  std::vector<double> distance;
  std::vector<double> psik;
  std::vector<double> electronDensity;
  std::vector<double> H2Density;
  std::vector<double> HIDensity;
};

class galAxions{
 private:
  
  double axionMassInEv; 
  double gagInGeV;
  double PagPDF;
  double IavPDF;
  double QavPDF;
  double UavPDF;
  double VavPDF;
  double PolDegPDF;
  double PosAnglePDF;
  
  std::string initFilename;
  std::string outputFilename;
  std::string galaxyFilename;
  
  std::ofstream outputStream;
  std::ofstream galaxyStream;
  
  MagneticField * magneticField; // = NULL;
  Ferriere * gas;
  
  los_struct los;
  
  std::vector<double> DeltaAgamma;
  std::vector<double> DeltaPl;
  std::vector<double> DeltaQED;
  std::vector<double> DeltaPar;
  std::vector<double> DeltaPerp;
  
 public:
  
  galAxions(const double& axionMassInEv_, 
	    const double& gagInGeV_,
	    const std::string& _initFilename)
    { 
      axionMassInEv = axionMassInEv_; 
      gagInGeV = gagInGeV_;
      initFilename = _initFilename;
      outputFilename = _initFilename; 
      outputFilename += ".txt";
      galaxyFilename = _initFilename;
      galaxyFilename += ".gal";
      //magnfield = NULL;
      
      outputStream.open(outputFilename.c_str(), std::ofstream::out);
      galaxyStream.open(galaxyFilename.c_str(), std::ofstream::out);
    }
  
  ~galAxions()
    {
      if (magneticField) delete magneticField;
      if (gas) delete gas;
      outputStream.close();
      galaxyStream.close();
    }
  
  double getAxionMassInEv(){ return axionMassInEv; }
  double getGagInGeV(){ return gagInGeV; }
  std::string getInitFilename(){ return initFilename; }
  
  void createGasDensity(void);  
  void createMagneticField(const MagneticFieldType& btype_,const int& bmode_=ASS);
  void createLos(const double&,const double&,const double&,const bool& printLos_=false);
  void printLos(const double&);
  std::vector<double> calculateProbability(const unsigned int&,const double&,const double&,const bool&,const bool&);
  std::vector<double> calculateProbability(const double&,const bool&,const bool&);
};

#endif