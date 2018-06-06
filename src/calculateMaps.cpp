#include <fstream>
#include <iostream>
#include <string>

#include "galaxions.h"
#include "healpixMap.h"

#define KEV 1e3
#define MEV 1e6
#define GEV 1e9
#define TEV 1e12

int main(int argc, char** argv)
{
//  std::cout<<"#Welcome to GalAxions Maps!"<<std::endl;
//
//  if (argc != 5) {
//    std::cerr<<"Usage: ./calculateMaps.out <output file name> mass [eV] coupling [GeV^-1] obs Energy [GeV]"<<std::endl;
//    std::cerr<<"Example: ./calculateMaps.out 1e-13 1e-11 1"<<std::endl;
//    return -1;
//  }
//
//  //const double mass = atof(argv[2])/1e-13;
//  //const double gag = atof(argv[3])/1e-11;
//
//  const double axionMassInEv = atof(argv[2]); // massInEv;
//  const double gagInGeV = atof(argv[3]); // couplingInGeV;
//  const double observedEnergy = atof(argv[4]);
//
//  std::ostringstream s;
//  s<<argv[1]<<"_"<<std::scientific<<std::setprecision(2)<<axionMassInEv<<"_"<<gagInGeV<<"_"<<observedEnergy;
//
//  const std::string initFilename = s.str();
//  const std::string asciiFilename = s.str()+".maps";
//
//  const int healpixResolution = 6/*7*/;
//
//  healpixMap * h = new healpixMap (healpixResolution,initFilename);
//
//  float * IavPDFarray = new float[ h->getMaxIter() ];
//  double * ldegarray = new double[ h->getMaxIter() ];
//  double * bdegarray = new double[ h->getMaxIter() ];
//
//#ifdef _OPENMP
//#pragma omp parallel for ordered schedule(dynamic) default(shared)
//#endif
//  for (unsigned long counter = 0; counter < h->getMaxIter(); counter++) {
//
//    galAxions * g = new galAxions (axionMassInEv,gagInGeV,initFilename);
//
//    g->createMagneticField(FARRAR);
//    g->createGasDensity();
//
//    double l = 0;
//    double b = 0;
//
//    pix2ang_ring(h->getNside(),counter,&b,&l);
//    b = M_PI/2.0-b;
//
//    const double bdeg = b/DegToRad;
//    const double ldeg = l/DegToRad;
//
//    g->createLos(ldeg,bdeg,30.0);
//
//    std::vector<double> result = g->calculateProbability(observedEnergy*GEV,false,false);
//
//    bdegarray[counter] = bdeg;
//    ldegarray[counter] = ldeg;
//    IavPDFarray[counter] = result[1];
//
//    //std::cout<<counter<<"\t"<<bdeg<<"\t"<<ldeg<<"\t"<<result[0]<<"\t"<<result[1]<<"\t"<<result[2]<<std::endl;
//
//    if (g) delete g;
//  }
//
//  std::ofstream outfile( asciiFilename.c_str() );
//  for ( unsigned long counter = 0; counter < h->getMaxIter(); counter++ ){
//    outfile<<counter<<"\t"<<ldegarray[counter]<<"\t"<<bdegarray[counter]<<"\t"<<IavPDFarray[counter]<<std::endl;
//  }
//  outfile.close();
//
//  long nside = h->getNside();
//  std::string filename = h->getProbabilityMapFilename();
//
//  write_healpix_map(IavPDFarray, nside, filename.c_str(), '0', "G");
//  // float *signal, long nside, char *filename, char nest, char *coordsys
//
//  if (h) delete h;
//  return 0;
}
