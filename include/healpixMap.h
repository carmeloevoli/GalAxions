#ifndef __HEALPIXMAP_H
#define __HEALPIXMAP_H

#include "chealpix.h"

class healpixMap{
 private:
  
  unsigned long healpixResolution = 0;
  unsigned long nside = 0;
  unsigned long npix = 0;
  
  double resolutionInDegree = 0;

  std::string probabilityMapFilename; // FITS file
  std::string pdMapFilename; // FITS file
  std::string paMapFilename; // FITS file

 public:
  
  healpixMap() { }
  
  healpixMap(const unsigned long& _hResolution,
	     const std::string& _initFilename)
    {
      healpixResolution = _hResolution;
      nside = pow(2.0, healpixResolution);
      npix = 12*nside*nside;
      resolutionInDegree = 4.0*M_PI/sqrt((double)npix);

      probabilityMapFilename = "!";
	probabilityMapFilename += _initFilename;
      probabilityMapFilename += ".fits";
      pdMapFilename = _initFilename;
      pdMapFilename += "_PD.fits"; 
      paMapFilename = _initFilename;
      paMapFilename += "_PA.fits";
    }
  
  virtual ~healpixMap() { }
  
  inline unsigned long getNpix() const { return npix; }
  inline unsigned long getNside() const { return nside; }
  inline unsigned long getMaxIter() { return npix; }
  inline double getResolution() const {return resolutionInDegree; }  

  std::string getProbabilityMapFilename() const { return probabilityMapFilename; }
};

#endif
