#ifndef BETA_HH
#define BETA_HH

//////////////////////////////////////////////////
// Beta.hh
// 
// Class for beta properties/ kinematics (instead of letting variables float around)
//////////////////////////////////////////////////

namespace urchin
{
    class Beta
    {
        public:
            double x0; //Guiding center position (x0,y0 -> rho)
            double y0;
            double rho;
            double mass;
            double fc;
            double magnetic_field;
            double omega_c;
            double q; //absolute value of charge
            double k;
            double energy_J; 
            double energy_eV; 
            double beta;
            double velocity;
            double gamma;
            double Rc; //cyclotron radius

            Beta(const double &aFc, const double &aMagneticField, const double &aX0, const double &aY0);
            Beta(const double &aFc, const double &aMagneticField, const double &aX0);
            
            double LarmorPower();
    };
} /* namespace urchin */
#endif
