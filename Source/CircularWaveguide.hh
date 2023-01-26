#ifndef CIRCULARWAVEGUIDE_HH
#define CIRCULARWAVEGUIDE_HH
#include <string>
#include <vector>
#include "Waveguide.hh"

//////////////////////////////////////////////////
// CircularWaveguide.hh
// Derived class for circular waveguides
// Encapsulates details for power in each mode (n,m,h)
// To avoid recalculating Bessel zeros, reads in 1000x1000 array from storage
//////////////////////////////////////////////////

class CircularWaveguide: public Waveguide
{
    public:
        CircularWaveguide(const double &aA, const unsigned &aNRoots);
        double TEModePower(const int& n, const int& m, const int& h, const Beta& beta) override;
        double TMModePower(const int& n, const int& m, const int& h, const Beta& beta) override;
        double TEkc(const int& n, const int& m) override;
        double TMkc(const int& n, const int& m) override;

    private:
        double a; //waveguide radius
        unsigned nRoots; //number of zeros to store
        std::vector<std::vector<double> > jn_zeros;
        std::vector<std::vector<double> > jnp_zeros;
        std::vector<std::vector<double> > ReadCSV(const std::string &filename);
};
#endif
