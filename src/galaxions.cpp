#include "galaxions.h"

void galAxions::createMagneticField(const MagneticFieldType& btype_, const int& bmode_) {
	if (btype_ == CONSTANT) {
		std::cout << "# init constant magnetic field!" << std::endl;
		magneticField = std::make_shared<ConstantField>(2.0);
	} else if (btype_ == FARRAR) {
		std::cout << "# init Farrar magnetic field!" << std::endl;
		magneticField = std::make_shared<FarrarField>();
	} else if (btype_ == PSHIRKOV) {
		std::cout << "# init Pshirkov magnetic field!" << std::endl;
		magneticField = std::make_shared<PshirkovField>(bmode_, 0.0, 8.5, 2.0);
	} else {
		std::cerr << "Magnetic Field model not found!\n";
	}
}

void galAxions::createGasDensity(const GasDensityType& gastype_) {
	if (gastype_ == FERRIERE) {
		std::cout << "# init Ferriere gas" << std::endl;
		gas = std::make_shared<Ferriere>();
	} else {
		std::cerr << "Gas Density model not found!\n";
	}

}

void galAxions::printLos(const double& rmax_) {
	galaxyStream << std::scientific << std::setprecision(5);
	for (size_t i = 0; i < los.nSteps; i++) {
		if (los.distance[i] < rmax_) {
			galaxyStream << los.distance[i] / kpc << "\t";
			galaxyStream << los.magneticFieldPerp[i] << "\t";
			galaxyStream << los.magneticFieldTotal[i] << "\t";
			galaxyStream << los.psik[i] << "\t";
			galaxyStream << los.electronDensity[i] << "\t";
			galaxyStream << std::endl;
		}
	}
}

void galAxions::createLos(const double& ldeg_, const double& bdeg_, const double& maxDistance_) {
	los.ldeg = ldeg_;
	los.bdeg = bdeg_;

	los.l = los.ldeg * DegToRad; // \phi in [0,2\pi]
	//los.b = (bdeg < 0) ? (90.0-bdeg)*DegToRad : -1; // \theta in [0,\pi]
	los.b = los.bdeg * DegToRad;

	double sinl = sin(los.l);
	double cosl = cos(los.l);

	double sinb = sin(los.b);
	double cosb = cos(los.b);

	double cosbcosl = cosb * cosl;
	double cosbsinl = cosb * sinl;

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

	while (fabs(xGalactoCentric) < x_max && fabs(yGalactoCentric) < y_max && fabs(zGalactoCentric) < z_max && distanceAlongLos < maxDistance_) {

		distanceAlongLos += step_size;

		if (sun_x < 0)
			xGalactoCentric = sun_x + distanceAlongLos * cosbcosl;
		else
			xGalactoCentric = sun_x - distanceAlongLos * cosbcosl;

		yGalactoCentric = distanceAlongLos * cosbsinl;
		zGalactoCentric = distanceAlongLos * sinb;

		Bperp.clear();
		Bperp = magneticField->GetBperp(xGalactoCentric, yGalactoCentric, zGalactoCentric); // [muG]
		Btotal = magneticField->GetB(xGalactoCentric, yGalactoCentric, zGalactoCentric); // [muG]

		los.magneticFieldPerp.push_back(normVector(Bperp));
		los.magneticFieldTotal.push_back(normVector(Btotal));

		double cosarg = 0;
		double testpsik = 0; // the angle between Btransverse and y axis

		if (!done) { // Fix once the reference direction

			done = true;

			ydir = Bperp;
			ydir[0] /= los.magneticFieldPerp.back();
			ydir[1] /= los.magneticFieldPerp.back();
			ydir[2] /= los.magneticFieldPerp.back();

			/*xdir[0] = (ydir[1]*dir[2]-dir[1]*ydir[2]);
			 xdir[1] = -(ydir[0]*dir[2]-dir[0]*ydir[2]);
			 xdir[2] = (ydir[0]*dir[1]-dir[0]*ydir[1]);*/
		}

		if (los.magneticFieldPerp.back() != 0) {
			cosarg = (ydir[0] * Bperp[0] + ydir[1] * Bperp[1] + ydir[2] * Bperp[2]) / los.magneticFieldPerp.back();
			if (cosarg >= 0.99999999999)
				testpsik = 0.0;
			else if (cosarg <= -0.99999999999)
				testpsik = PI;
			else
				testpsik = acos(cosarg);
		} else
			testpsik = 0.0;

		los.psik.push_back(testpsik);
		los.distance.push_back(distanceAlongLos);
		los.electronDensity.push_back(gas->get(xGalactoCentric, yGalactoCentric, zGalactoCentric));
	}

	los.nSteps = los.distance.size();

	std::reverse(los.distance.begin(), los.distance.end());
	std::reverse(los.magneticFieldPerp.begin(), los.magneticFieldPerp.end());
	std::reverse(los.magneticFieldTotal.begin(), los.magneticFieldTotal.end());
	std::reverse(los.psik.begin(), los.psik.end());
	std::reverse(los.electronDensity.begin(), los.electronDensity.end());
}
