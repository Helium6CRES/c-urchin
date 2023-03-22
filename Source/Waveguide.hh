#ifndef WAVEGUIDE_HH
#define WAVEGUIDE_HH
#include <fstream>
#include <string>

#include "Beta.hh"

namespace urchin
{

    class Waveguide
    {
        public:
            virtual ~Waveguide() = default;

            virtual double TEModePower(const int &n, const int &m, const int &h, const Beta &beta) = 0;
            virtual double TMModePower(const int &n, const int &m, const int &h, const Beta &beta) = 0;

            virtual bool HitsWall(const Beta &beta) = 0;

            double ModePower(const int &n, const int &m, const int &h, const Beta &beta, const bool &bTE);

            virtual double TEkc(const int &n, const int &m) = 0;
            virtual double TMkc(const int &n, const int &m) = 0;
            double kc(const int &n, const int &m, const bool &bTE);
            
            double TotalPower(const unsigned &N, const unsigned &M, const unsigned &H, const double &fTolerance, const Beta &beta, const bool &bTE);
            double TotalPowerHarmonic(const unsigned &N, const unsigned &M, const unsigned &h, const Beta &beta, const bool &bTE);


            void OpenCSV(const std::string &filename);
            void WriteCSV(const Beta &beta, const double &totalPower, const double &te11Power);
            void CloseCSV();

            std::ofstream file;
    };

} /* namespace urchin */
#endif
