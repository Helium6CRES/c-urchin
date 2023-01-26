#include <boost/math/special_functions/bessel.hpp> 
#include <boost/math/special_functions/bessel_prime.hpp> 
#include <fstream>
#include <sstream>
#include <iostream>

#include "CircularWaveguide.hh"
#include "Constants.hh"
#include "Param.hh"

CircularWaveguide::CircularWaveguide(const double &aA, const unsigned &aNRoots):
        a(aA),
        nRoots(aNRoots), 
        jn_zeros( nRoots,std::vector<double> (nRoots)), 
        jnp_zeros( nRoots,std::vector<double> (nRoots))
{
    jn_zeros = ReadCSV(DATA_DIR + std::string("/jns.csv"));
    jnp_zeros = ReadCSV(DATA_DIR + std::string("/jnps.csv"));

    return;
}

std::vector<std::vector<double> > CircularWaveguide::ReadCSV(const std::string &filename)
{
    std::ifstream inputFile(filename);
    std::istringstream in;
    if(!inputFile.is_open()) throw std::runtime_error("Could not open bessel zeros file!");

    std::vector<std::vector<double> > zerosVector;
    std::vector<double> vTmp;
    std::string line;

    while(std::getline(inputFile,line))
    {
        in.clear();
        in.str(line);
        vTmp.clear();
        for(double value; in >> value; in.ignore())
        {
            vTmp.push_back(value);
        }

        zerosVector.push_back(vTmp);
    }

    inputFile.close();

    return zerosVector;
}

double CircularWaveguide::TEkc(const int& n, const int& m)
{
    return jnp_zeros[n][m] / a;
}

double CircularWaveguide::TMkc(const int& n, const int& m)
{
    return jn_zeros[n][m] / a;
}
double CircularWaveguide::TEModePower(const int& n, const int& m, const int& h, const Beta& beta)
{
    if(!m) return 0;
    const double p_prime = jnp_zeros[n][m];
    const double kc = TEkc(n,m);
    double output = pow(boost::math::cyl_bessel_j(n+h, kc*beta.rho),2) + pow(boost::math::cyl_bessel_j(abs(n-h), kc*beta.rho),2);
    output *= pow(boost::math::cyl_bessel_j_prime(h, kc*beta.Rc),2);

    output *= pow(beta.q*beta.velocity * h * beta.omega_c * Constants::MuNull() / kc,2);
    const double wg_beta = sqrt(pow(h*beta.k,2)-pow(kc,2));
    const double P_N = (Constants::Pi() * h * beta.omega_c * Constants::MuNull() * wg_beta ) / (2*pow(kc,4))*(pow(p_prime,2)-n*n) * pow(boost::math::cyl_bessel_j(n,p_prime),2);

    return output / (2*P_N);
}

double CircularWaveguide::TMModePower(const int& n, const int& m, const int& h, const Beta& beta)
{
    if(!m) return 0;
    const double p_nm = jn_zeros[n][m];
    const double kc = TMkc(n,m);
    double output = pow(boost::math::cyl_bessel_j(n+h, kc*beta.rho),2) + pow(boost::math::cyl_bessel_j(abs(n-h), kc*beta.rho),2);
    output *= pow(boost::math::cyl_bessel_j(h, kc*beta.Rc),2);

    const double wg_beta = sqrt(pow(h*beta.k,2)-pow(kc,2));
    output *= pow(beta.q*beta.velocity * h * wg_beta / (kc*kc*beta.Rc) ,2);
    const double P_N = (Constants::Pi() * h * beta.omega_c * Constants::EpsNull() * wg_beta ) / (2*pow(kc,4)) * pow(p_nm,2) * pow(boost::math::cyl_bessel_j_prime(n,p_nm),2);

    return output / (2*P_N);
}
