#include <iostream>
#include <string>

#include "galaxions.h"

int main(int argc, char** argv)
{
  std::cout<<"#Welcome to GalAxions!"<<std::endl; 
  
  /*if (argc != 4) {
    cerr << "Usage: ./pagamma <output file name> mass [eV] coupling [GeV^-1]" << endl;
    return -1;
    }*/
  
  const double mass = 1e-13; // massInEv;
  const double gag  = 1e-11; // couplingInGeV;
  const std::string initFilename = "temp";
  
  galAxions * g = new galAxions (mass,gag,initFilename);
  
  g->createMagneticField(FARRAR);
  
  g->createGasDensity();

  g->createLos(279.6,-32.1); // SN1987A
      
  //const double mass = atof(argv[2])/1e-13;
  //const double gag = atof(argv[3])/1e-11;
  
  if (g) delete g; 
  return 0;
}
