#include <iostream>
#include <string>

#include "galaxions.h"

#define KEV 1e3
#define MEV 1e6
#define GEV 1e9
#define TEV 1e12

int main(int argc, char** argv)
{
  std::cout<<"#Welcome to GalAxions!"<<std::endl; 
  
  if (argc != 4) {
    std::cerr<<"Usage: ./galAxions.exe <output file name> mass [eV] coupling [GeV^-1]"<<std::endl;
    return -1;
  }
  
  //const double mass = atof(argv[2])/1e-13;
  //const double gag = atof(argv[3])/1e-11;
  
  const double axionMassInEv = atof(argv[2]); // massInEv;
  const double gagInGeV = atof(argv[3]); // couplingInGeV;
  
  const std::string initFilename = argv[1];
  
  galAxions * g = new galAxions (axionMassInEv,gagInGeV,initFilename);
  
  g->createMagneticField(FARRAR);
  //g->createMagneticField(PSHIRKOV,ASS);
  //g->createMagneticField(CONSTANT);
  
  g->createGasDensity();
  
  //g->createLos(279.6,-32.1,30); // SN1987A
  //g->createLos(0.01,0.01,10); // GC
  g->createLos(338.32,-0.02,8.6); // HESS J1640-465
  
  //g->calculateProbability(100,1.0*KEV,1e2*KEV,false,true);
  //g->calculateProbability(2000,1.0*MEV,100.0*MEV,false,true); // Payez Model
  //g->calculateProbability(2000,1.0*MEV,300.0*MEV,false,true); // Payez Model GC
  g->calculateProbability(2000,1e-4*TEV,1e2*TEV,false,true); // HESS J1640-465

  if (g) delete g; 
  return 0;
}
