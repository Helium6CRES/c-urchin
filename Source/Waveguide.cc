#include <iomanip>
#include <iostream>
#include <math.h>
#include <boost/filesystem.hpp>

#include "Waveguide.hh"

namespace urchin
{

    double Waveguide::ModePower(const int &n, const int &m, const int &h, const Beta &beta, const bool &bTE) 
    {
        return bTE? TEModePower(n,m,h,beta): TMModePower(n,m,h,beta);
    }

    double Waveguide::kc(const int &n, const int &m, const bool &bTE)
    {
        return bTE? TEkc(n,m): TMkc(n,m);
    }

    double Waveguide::TotalPower(const unsigned &N, const unsigned &M, const unsigned &H, const double &fTolerance, const Beta& beta, const bool &bTE)
    {
        double totalSum = 0.;
        const bool bfTol = bool(fTolerance);

        double hSum, tmpModePower;

        for(int n=0; n < N; ++n)
            for(int m=0; m < M; ++m)
            {
                int h_min = int(ceil(kc(n,m,bTE)/beta.k));
                hSum = 0.;

                for(int h=h_min; h < H; ++h)
                {
                    tmpModePower = ModePower(n, m, h, beta,bTE);
                    hSum += tmpModePower;
                    if(bfTol && tmpModePower / hSum < fTolerance) break;
                }

                totalSum +=hSum;
            }

        return totalSum;
    }

    void Waveguide::OpenCSV(const std::string &filename)
    {
        const std::string dlmr = ","; //delimeter
        bool bFileExists = bool(boost::filesystem::exists(filename));
        file.open(filename, std::ios::app);
        if(!bFileExists)
            file <<"frequency"<<dlmr<<"magnetic_field"<<dlmr<<"rho"<<dlmr<<"power"<<"\n";

         file << std::setprecision(16);
    }

    void Waveguide::WriteCSV(const Beta &beta, const double &power)
    {
        const std::string dlmr = ","; //delimeter
        file << beta.fc<<dlmr<<beta.magnetic_field<<dlmr<<beta.rho<<dlmr<<power<<"\n";
    }

    void Waveguide::CloseCSV()
    {
        file.close();
    }

} /* namespace urchin */
