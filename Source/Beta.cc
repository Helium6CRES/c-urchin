#include <math.h>

#include "Constants.hh"
#include "Beta.hh"

namespace urchin
{

    Beta::Beta(const double &aFc, const double &aMagneticField, const double &aX0, const double &aY0):
            fc(aFc),
            magnetic_field(aMagneticField),
            x0(aX0),
            y0(aY0)
    {
        rho = sqrt(x0 * x0 + y0 * y0);
        const double c = Constants::C();
        q = Constants::Q();
        mass = Constants::M_el_kg();
        omega_c = 2. * Constants::Pi() * fc;
        gamma = q * magnetic_field / (mass * omega_c);
        beta = sqrt(1-pow(1./gamma,2));
        velocity = beta * c;
        energy_J = (gamma-1)* mass * c*c;
        energy_eV = energy_J / q;
        k = omega_c / c;
        Rc = velocity / omega_c;
    }

    //Delegating constructor if y0 not passed
    Beta::Beta(const double &aFc, const double &aMagneticField, const double &aX0): Beta(aFc,aMagneticField, aX0, 0) {};

    double Beta::LarmorPower()
    {
        return q*q*Constants::C()*pow(beta*gamma,4) / ( 6. * Constants::Pi() * Constants::EpsNull() * Rc * Rc);
    }

} /* namespace urchin */
