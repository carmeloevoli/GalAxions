#ifndef __GALAXIONS_H
#define __GALAXIONS_H

#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <memory>
#include <fstream>
#include <string>
#include <complex>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gas.h"
#include "Ferriere07.h"
#include "magneticfield.h"
#include "Farrar.h"
#include "Pshirkov.h"
#include "JF12.h"
#include "mymatrix.h"
#include "solvers.h"

enum MagneticFieldType {
	CONSTANT, FARRAR, PSHIRKOV
};

enum GasDensityType {
	FERRIERE, CORDES, YMW
};

struct LOS {
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
};

class galAxions {
private:

	double axionMass;
	double gag;

	std::string outputFilename;
	std::string losFilename;
	std::ofstream outputStream;
	std::ofstream galaxyStream;

	std::shared_ptr<MagneticField> magneticField;
	std::shared_ptr<ElectronModel> gas;

	LOS los;

	/*double PagPDF;
	 double IavPDF;
	 double QavPDF;
	 double UavPDF;
	 double VavPDF;
	 double PolDegPDF;
	 double PosAnglePDF;
	 */
	/*std::vector<double> DeltaAgamma;
	 std::vector<double> DeltaPl;
	 std::vector<double> DeltaQED;
	 std::vector<double> DeltaPar;
	 std::vector<double> DeltaPerp;*/

public:

	galAxions(const double& axionMass_, const double& gag_, const std::string& initFilename_) :
			axionMass(axionMass_), gag(gag_) {
		outputFilename = initFilename_ + ".txt";
		losFilename = initFilename_ + ".los";
		outputStream.open(outputFilename.c_str(), std::ofstream::out);
		galaxyStream.open(losFilename.c_str(), std::ofstream::out);
	}

	~galAxions() {
		outputStream.close();
		galaxyStream.close();
	}

	void createGasDensity(const GasDensityType& gastype_);

	void createMagneticField(const MagneticFieldType& btype_, const int& bmode_ = ASS);

	void createLos(const double& ldeg_, const double& bdeg_, const double& maxDistance_ = 30 * kpc);

	void printLos(const double&);

	void reverseLos();

	std::vector<double> calculateProbability(const unsigned int&, const double&, const double&, const bool&, const bool&);

	std::vector<double> calculateProbability(const double&, const bool&, const bool&);
};

#endif
