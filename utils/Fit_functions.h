// Fitting functions

#define PionPdgMass        0.139570
#define ProtonPdgMass      0.938272
#define KaonPdgMass        0.493677  //charged kaons, not K0s
#define RhoPdgMass         0.770
#define RhoLifeTime        0.148
#define OmegaPdgMass       0.782
#define OmegaLifeTime      0.017

#include <TMath.h>
#include <complex>
#include <cmath>

std::complex<double> pole_term(double q, double M, double Gamma) {
    /**
    * @brief Characteristic pole term of a particle's propagator.
    * 
    * @param q : momentum.
    * @param M : Particle's mass.
    * @param Gamma : Particle's decay rate.
    * @return std::complex<double> 
    */

    double numerator = sqrt(q * M * Gamma);
    std::complex<double> denominator(q * q - M * M, + Gamma * M);

    return numerator/denominator;
}

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

    std::complex<double> rho_term = A_rho*pole_term(M, M_rho, Gamma_rho);
    complex<double> exp_w(TMath::Cos(phi_omega), TMath::Sin(phi_omega));
    std::complex<double> omega_term = C_omega*pole_term(M, M_omega, Gamma_omega)*exp_w;

    return std::norm(rho_term+omega_term+B_pp);
}

// Breit-Wigner Function + Soding Term.
//

double Pole_term_squared(const double *x, const double *par) {
    double A =  par[0];
    double B = par[1];
    double M_rho = par[2];
    double Gamma_rho = par[3];

    double M = x[0];
    return std::norm(A*pole_term(M, M_rho, Gamma_rho));
}

double Interference_term(const double *x, const double *par) {
    double A =  par[0];
    double B = par[1];
    double M_rho = par[2];
    double Gamma_rho = par[3];

    double M = x[0];
    return 2*B*(A*pole_term(M, M_rho, Gamma_rho)).real();
}

double SodingEqn(const double *x, const double *par) {
    double A =  par[0];
    double B = par[1];
    double M_rho = par[2];
    double Gamma_rho = par[3];
    double M = x[0];
    return std::norm(A*pole_term(M, M_rho, Gamma_rho) + B);
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
    return 1.0 - par[0] * TMath::Cos(2.0 * x[0]);
}