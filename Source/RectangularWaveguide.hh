#ifndef RECTANGULARWAVEGUIDE_HH
#define RECTANGULARWAVEGUIDE_HH
#include "Waveguide.hh"

//////////////////////////////////////////////////
// RectangularWaveguide.hh
// Derived class for rectangular waveguides
// Encapsulates details for power in each mode (n,m,h)
//////////////////////////////////////////////////

class RectangularWaveguide: public Waveguide
{
    public:
        RectangularWaveguide(const double &aA, const double &aB);
        double TEModePower(const int& n, const int& m, const int& h, const Beta& beta) override;
        double TMModePower(const int& n, const int& m, const int& h, const Beta& beta) override;
        double TEkc(const int& n, const int& m) override;
        double TMkc(const int& n, const int& m) override;

    private:
        double a;  //waveguide dimensions
        double b; 

        double kx(const int& m);
        double ky(const int& n);
};
#endif
