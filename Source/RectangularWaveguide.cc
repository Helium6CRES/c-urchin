#include <boost/math/special_functions/bessel.hpp> 
#include <boost/math/special_functions/bessel_prime.hpp> 

#include "Constants.hh"
#include "RectangularWaveguide.hh"

RectangularWaveguide::RectangularWaveguide(const double &aA, const double &aB): a(aA), b(aB) {}

double RectangularWaveguide::kx(const int& m)
{
    return m * Constants::Pi() / a;
}

double RectangularWaveguide::ky(const int& n)
{
    return n * Constants::Pi() / b;
}

double RectangularWaveguide::TEkc(const int& n, const int& m)
{
    const double tKx = kx(m);
    const double tKy = ky(n);
    return sqrt( tKx * tKx + tKy * tKy);
}

double RectangularWaveguide::TMkc(const int& n, const int& m)
{
    return TEkc(n,m);
}

double RectangularWaveguide::TEModePower(const int& n, const int& m, const int& h, const Beta& beta)
{
    if(!m && !n) return 0;
    const double kc = TEkc(n,m);
    const double delta = atan2(-n/b,m/a);
    const double tKx = kx(m);
    const double tKy = ky(n);
    
    const double phix0 = kx(m) * beta.x0;
    const double phiy0 = ky(n) * beta.y0;

    double output;

    if(h%2) //odd
    {
        output = pow(sin(phix0) * cos(phiy0) *cos(h * delta),2) + pow(cos(phix0) * sin(phiy0) *sin(h * delta),2);
    }
    else
    {
        output = pow(cos(phix0) * cos(phiy0) *cos(h * delta),2) + pow(sin(phix0) * sin(phiy0) *sin(h * delta),2);
    }

    output *= pow(beta.q * beta.velocity * h * beta.omega_c * Constants::MuNull() * boost::math::cyl_bessel_j_prime(h,kc*beta.Rc) / kc,2);
    const double wg_beta = sqrt(pow(h*beta.k,2)-pow(kc,2));
    double P_N = h * beta.omega_c * Constants::MuNull() * wg_beta * a * b / (8 * pow(kc,4));
    const double dn = double(!n);
    const double dm = double(!m);
    P_N *= (tKx * tKx * (1-dm) * (1+dn)   + tKy * tKy * (1+dm) * (1-dn));
    return output / P_N;
}

double RectangularWaveguide::TMModePower(const int& n, const int& m, const int& h, const Beta& beta)
{
    if(!m && !n) return 0;
    const double kc = TEkc(n,m);
    const double delta = atan2(-n/b,m/a);
    const double tKx = kx(m);
    const double tKy = ky(n);
    
    const double phix0 = kx(m) * beta.x0;
    const double phiy0 = ky(n) * beta.y0;

    double output;

    if(h%2) //odd
    {
        output = pow(cos(phix0) * sin(phiy0) *cos(h * delta),2) + pow(sin(phix0) * cos(phiy0) *sin(h * delta),2);
    }
    else
    {
        output = pow(sin(phix0) * sin(phiy0) *cos(h * delta),2) + pow(cos(phix0) * cos(phiy0) *sin(h * delta),2);
    }

    const double wg_beta = sqrt(pow(h*beta.k,2)-pow(kc,2));
    output *= pow(beta.q * beta.velocity * wg_beta* h * boost::math::cyl_bessel_j(h,kc*beta.Rc) /(kc * kc * beta.Rc) ,2);
    double P_N = h * beta.omega_c * Constants::EpsNull() * wg_beta * a * b / (8 * pow(kc,4));
    const double dn = double(!n);
    const double dm = double(!m);
    P_N *= (tKx * tKx * (1+dm) * (1-dn)   + tKy * tKy * (1-dm) * (1+dn));
    return output / P_N;
}

