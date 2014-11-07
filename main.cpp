#include <iostream>
#include <string>

#include "galaxions.h"

#define KEV 1e3
#define MEV 1e6
#define GEV 1e9

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
  
  g->createLos(279.6,-32.1); // SN1987A
  
  g->calculateProbability(100,1.0*KEV,1e2*KEV,false,true);
  //g->calculateProbability(10000,1,300.0*MEV,false); Payez Model

  if (g) delete g; 
  return 0;
}
