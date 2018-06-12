#include "galconda.h"

#define ALMOST_ONE 0.9999999999

void GALCONDA::createMagneticField(const MagneticFieldType& btype_, const int& bmode_) {
    if (btype_ == CONSTANT) {
        std::cout << "# init constant magnetic field!" << std::endl;
        magnetic_field = std::make_shared<ConstantField>(2.0);
    }
    else if (btype_ == FARRAR) {
        std::cout << "# init Farrar magnetic field!" << std::endl;
        magnetic_field = std::make_shared<JF12Field>();
    }
    else if (btype_ == PSHIRKOV) {
        std::cout << "# init Pshirkov magnetic field!" << std::endl;
        magnetic_field = std::make_shared<PshirkovField>(bmode_, 0.0, 8.5, 2.0);
    }
    else {
        std::cerr << "Magnetic Field model not found!\n";
    }
}

void GALCONDA::createGasDensity(const GasDensityType& gastype_) {
    if (gastype_ == FERRIERE) {
        std::cout << "# init Ferriere gas\n";
        gas = std::make_shared<Ferriere07>();
    }
    else if (gastype_ == CORDES) {
        std::cout << "# init Cordes gas\n";
        gas = std::make_shared<Cordes91>();
    }
    else if (gastype_ == YMW) {
        std::cout << "# init YMW gas\n";
        gas = std::make_shared<YMW16>();
    }
    else {
        std::cerr << "Gas Density model not found!\n";
    }

}

void GALCONDA::printLos(const double& rmax_) {
    std::cout << "# write los in " << los_filename << "\n";
    los_ss << "# d [kpc] - B perp [muG] - B tot [muG] - psi_k - n_e [cm-3]\n";
    los_ss << std::scientific << std::setprecision(5);
    for (auto domain : los.domains) {
        los_ss << domain.distance / kpc << "\t";
        los_ss << domain.magnetic_field_perp / muG << "\t";
        los_ss << domain.magnetic_field_total / muG << "\t";
        los_ss << domain.psik << "\t";
        los_ss << domain.electron_density * cm3 << "\t";
        los_ss << std::endl;
    }
}

double get_psik(const std::vector<double>& ref_direction, const std::vector<double>& B_perp) {
    double psik = 0; // the angle between B_transverse and y axis
    auto B_perp_norm = normVector(B_perp);
    if (B_perp_norm != 0) {
        double cosarg = (ref_direction[0] * B_perp[0] + ref_direction[1] * B_perp[1] + ref_direction[2] * B_perp[2]);
        cosarg /= B_perp_norm;
        if (cosarg >= ALMOST_ONE)
            psik = 0.0;
        else if (cosarg <= -ALMOST_ONE)
            psik = M_PI;
        else
            psik = std::acos(cosarg);
    }
    else {
        psik = 0.0;
    }
    return psik;
}

void GALCONDA::createLos(const double& ldeg_, const double& bdeg_, const double& maxDistance_) {
    std::cout << "# create los along l : " << ldeg_ << " and b : " << bdeg_ << "\n";

    los = LOS(ldeg_, bdeg_);

    double distance_along_los = 0.0;
    double x_galaxy = 0.0;
    double y_galaxy = 0.0;
    double z_galaxy = 0.0;

    std::vector<double> ref_direction;

    bool done = false;

    while (fabs(x_galaxy) < x_max && fabs(y_galaxy) < y_max && fabs(z_galaxy) < z_max
            && distance_along_los < maxDistance_) {

        distance_along_los += step_size;

        if (sun_x < 0)
            x_galaxy = sun_x + distance_along_los * los.cos_b_cos_l;
        else
            x_galaxy = sun_x - distance_along_los * los.cos_b_cos_l;
        y_galaxy = distance_along_los * los.cos_b_sin_l;
        z_galaxy = distance_along_los * los.sin_b;

        auto B_perp = magnetic_field->GetBperp(x_galaxy, y_galaxy, z_galaxy);
        auto B_perp_norm = normVector(B_perp);
        auto B_total = magnetic_field->GetB(x_galaxy, y_galaxy, z_galaxy);
        auto B_total_norm = normVector(B_total);
        auto n_e = gas->get(x_galaxy, y_galaxy, z_galaxy);

        if (!done) { // Fix once the reference direction
            done = true;
            ref_direction = B_perp;
            ref_direction[0] /= B_perp_norm;
            ref_direction[1] /= B_perp_norm;
            ref_direction[2] /= B_perp_norm;
        }

        domain d = {B_perp_norm, B_total_norm, distance_along_los, get_psik(ref_direction, B_perp), n_e, 0, 0};

        los.domains.push_back(d);
    }

    std::reverse(los.domains.begin(), los.domains.end());
}

void GALCONDA::calculateProbability(const size_t& nEnergy_, const double& Emin_, const double& Emax_,
        const bool& do_damping_, const bool& do_output_) {

    time_t timeBegin, timeEnd;

    const int nDomains = los.domains.size();

    time(&timeBegin);

#ifdef _OPENMP
#pragma omp parallel for ordered schedule(dynamic) default(shared)
#endif
    for (size_t iE = 0; iE < nEnergy_; iE++) {

        const double Energy = (nEnergy_ == 1) ? Emin_ : pow(10,
                log10(Emin_) + double(iE) / double(nEnergy_ - 1) * log10(Emax_ / Emin_));

        std::vector<double> Delta_agamma;
        std::vector<double> Delta_pl;
        std::vector<double> Delta_QED;
        std::vector<double> Delta_par;
        std::vector<double> Delta_perp;

        const double Xsec_H2 = 0;
        const double Xsec_HI = 0;

        const double gag_norm = gag / (5e-11 / GeV);
        const double axion_mass_norm = axionMass / (1e-8 * eV);
        const double Energy_norm = Energy / TeV;
        const double Delta_a_kpc = -7.8e-3 * pow2(axion_mass_norm) / Energy_norm;

        for (auto domain : los.domains) {
            double B_T_norm = domain.magnetic_field_perp / muG;
            double n_e_norm = domain.electron_density / (1e-3 / cm3);
            double Delta_agamma_kpc = 7.6e-2 * gag_norm * B_T_norm; // Eq.4 in Horns+12
            double Delta_pl_kpc = -1.1e-10 / Energy_norm * n_e_norm;
            double Delta_QED_kpc = 4.1e-6 / Energy_norm * pow2(B_T_norm);

            Delta_agamma.push_back(Delta_agamma_kpc / kpc);
            Delta_pl.push_back(Delta_pl_kpc / kpc);
            Delta_QED.push_back(Delta_QED_kpc / kpc);
            Delta_par.push_back(Delta_pl.back() + 3.5 * Delta_QED.back()); // Eq.3.10 Bassan+10
            Delta_perp.push_back(Delta_pl.back() + 2.0 * Delta_QED.back()); // Eq.3.11 Bassan+10
        }

        MyMatrix Tk_total, rho, rho_old;

        rho_old(0, 0) = std::complex<double>(0.0, 0.0);
        rho_old(1, 1) = std::complex<double>(0.0, 0.0); // Photon unpolarized. Initial condition.
        rho_old(2, 2) = std::complex<double>(1.0, 0.0); // Full ALPs beam. Initial condition.

        Solver(nDomains, Delta_agamma, Delta_pl, Delta_QED, Delta_par, Delta_perp, Delta_a_kpc / kpc, los.domains,
                Xsec_H2, Xsec_HI, rho_old, rho, Tk_total);

        std::cout << "... energy : " << Energy / eV << std::endl;

        std::cout << "... T_k : " << std::endl;

        Tk_total.Print();

        double Pag_PDF = real(rho(2, 2));
        double Iav_PDF = real(rho(0, 0) + rho(1, 1));
        double Qav_PDF = real(rho(0, 0) - rho(1, 1));
        double Uav_PDF = real(rho(0, 1) + rho(1, 0));
        double Vav_PDF = imag(rho(1, 0) - rho(0, 1));

        double PolDeg_PDF =
                std::sqrt(pow2(Qav_PDF) + pow2(Uav_PDF) /*+ pow(VavPDF, 2)*/) / Iav_PDF; // Eq. 3.44 Bassan+10
        double PosAngle_PDF = 0.5 * atan2(Uav_PDF, Qav_PDF); // TODO convert in deg?

        if (Iav_PDF < 0) {
            std::cout << "Warning: Negative Intensity" << std::endl;
            rho.Print();
        }

        if (do_output_) {
            output_ss << std::scientific << Energy / eV << "\t" << Iav_PDF << "\t" << Pag_PDF << "\t";
            output_ss << real(rho(0, 0)) << "\t" << real(rho(1, 1)) << "\t" << real(rho(2, 2)) << "\t";
            output_ss << imag(rho(0, 0)) << "\t" << imag(rho(1, 1)) << "\t" << imag(rho(2, 2)) << "\t";
            output_ss << std::endl;
        }
    } //close energy loop

    time(&timeEnd);

    std::cout << "Ended in " << difftime(timeEnd, timeBegin) << " seconds." << std::endl;
}
