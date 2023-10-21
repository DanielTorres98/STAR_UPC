// Fitting functions

#include <TMath.h>
#include <complex>
#include <cmath>

double BW_rho_omega_photoproduction(const double *x, const double *par) {
    // Breit-Wigner function including rho and omega photoproduction
    //
    double M_rho = par[0];
    double Gamma_rho = par[1];
    double A_rho = par[2];
    double B_pp = par[3];
    double C_omega = par[4];
    double M_omega = par[5];
    double Gamma_omega = par[6];
    double phi_omega = par[7];
    double M = x[0];


    double numerator_rho = A_rho * sqrt(M * M_rho * Gamma_rho);
    std::complex<double> denominator_rho(M * M - M_rho * M_rho, + Gamma_rho * M_rho);
    std::complex<double> rho_term  = numerator_rho/denominator_rho;


    std::complex<double> numerator_omega(C_omega * sqrt(M * M_omega * Gamma_omega)*
                                        TMath::Cos(phi_omega), TMath::Sin(phi_omega));
    std::complex<double> denominator_omega(M * M - M_omega * M_omega, + Gamma_omega * M_omega);
    std::complex<double> omega_term  = numerator_omega/denominator_omega;

    double magnitude =  std::abs(rho_term + B_pp + omega_term);

    return magnitude*magnitude;
}


double SodingEqn(const double *x, const double *par) {
    double A =  par[0];
    double B = par[1];
    double M_rho = par[2];
    double Gamma_rho = par[3];

    double numerator = A * sqrt(x[0] * M_rho * Gamma_rho);
    std::complex<double> denominator(x[0] * x[0] - M_rho * M_rho, + Gamma_rho * M_rho);
    std::complex<double> denominatorCC(x[0] * x[0] - M_rho * M_rho, - Gamma_rho * M_rho);
    
    // Calculate the magnitude of the complex result
    double magnitude = std::abs((numerator / denominator + B)*(numerator / denominatorCC + B));
    
    // Return the magnitude as a double
    return magnitude;
}

double SodingEqn_addition_term(const double *x, const double *par) {
    double A =  par[0];
    double B = par[1];
    double M_rho = par[2];
    double Gamma_rho = par[3];

    double numerator = A * sqrt(x[0] * M_rho * Gamma_rho);
    std::complex<double> denominator(x[0] * x[0] - M_rho * M_rho, + Gamma_rho * M_rho);
    std::complex<double> denominatorCC(x[0] * x[0] - M_rho * M_rho, - Gamma_rho * M_rho);

    std::complex<double> Z(numerator/denominator);
    std::complex<double> ZCC(numerator/denominatorCC);

    return std::abs(Z*ZCC);
}

double SodingEqn_interference_term(const double *x, const double *par) {
    double A =  par[0];
    double B = par[1];
    double M_rho = par[2];
    double Gamma_rho = par[3];

    double numerator = A * sqrt(x[0] * M_rho * Gamma_rho);
    std::complex<double> denominator(x[0] * x[0] - M_rho * M_rho, + Gamma_rho * M_rho);
    std::complex<double> denominatorCC(x[0] * x[0] - M_rho * M_rho, - Gamma_rho * M_rho);

    std::complex<double> Z(numerator/denominator);
    std::complex<double> ZCC(numerator/denominatorCC);
    return B*(Z+ZCC).real();
}

double BW_Soding(const double *x, const double *par) {
    // Breit Wigner + Soding interference term.
    //

    double rhoMass = RhoPdgMass;
    double rhoLifetime = RhoLifeTime;
    
    double denominator = ((par[2]*par[2]-x[0]*x[0])*(par[2]*par[2]-x[0]*x[0])+
                         par[2]*par[2]*par[3]*par[3]);
    double BW = par[0]*x[0]*par[2]*par[3] / denominator;
    double Soding_Term = par[1]*(par[2]*par[2] - x[0]*x[0]) / denominator;

    return BW + Soding_Term;
}

double Ross_Stodolsky(const double *x, const double *par) {
    // Breit Wigner + Soding interference term.
    //

    double rhoMass = RhoPdgMass;
    double rhoLifetime = RhoLifeTime;
    
    double denominator = ((par[2]*par[2]-x[0]*x[0])*(par[2]*par[2]-x[0]*x[0])+
                         par[2]*par[2]*par[3]*par[3]);
    double BW = par[0]*x[0]*par[2]*par[3] / denominator;

    return BW*TMath::Power((rhoMass/x[0]), par[1]);
}

Double_t Interference_phi_fit(Double_t* x, Double_t* par) {
    return 1.0 - par[0] * cos(2.0 * x[0]);
}