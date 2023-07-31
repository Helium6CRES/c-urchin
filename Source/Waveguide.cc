#include <iomanip>
#include <iostream>
#include <math.h>
#include <boost/filesystem.hpp>

#include "Waveguide.hh"

namespace urchin
{

    double Waveguide::ModePower(const int &n, const int &m, const int &h, const Beta &beta, const bool &bTE) 
    {
        if(HitsWall(beta))
            return 0.;

        return bTE? TEModePower(n,m,h,beta): TMModePower(n,m,h,beta);
    }

    double Waveguide::kc(const int &n, const int &m, const bool &bTE)
    {
        return bTE? TEkc(n,m): TMkc(n,m);
    }

    double Waveguide::TotalPowerHarmonic(const unsigned &N, const unsigned &M, const unsigned &h, const Beta &beta, const bool &bTE)
    {
        if(HitsWall(beta))
            return 0.;

        double totalSum, tmpModePower;
        totalSum = 0;

        for(int n=0; n < N; ++n)
        {
            if( h <= kc(n,0,bTE)/beta.k) break;

            for(int m=0; m < M; ++m)
            {
                if( h <= kc(n,m,bTE)/beta.k) break;

                tmpModePower = ModePower(n, m, h, beta,bTE);
                totalSum +=tmpModePower;
            }
        }

        return totalSum;
    }

    double Waveguide::TotalPower(const unsigned &N, const unsigned &M, const unsigned &H, const double &fTolerance, const Beta& beta, const bool &bTE)
    {
        if(HitsWall(beta))
            return 0.;

        double totalSum = 0.;
        const bool bfTol = bool(fTolerance);

        double tmpModePower;
        double mSum, nSum;
        const int dH = 10;
        int h;

        for(int dh=0; dh < dH; ++dh)
        {
            nSum = 0.;
            for(int n=0; n < N; ++n)
            {
                mSum = 0.;
                for(int m=0; m < M; ++m)
                {
                    h = int(ceil(kc(n,m,bTE)/beta.k)) + dh;
                    if(h >= H) break;
                    tmpModePower = ModePower(n, m, h, beta,bTE);
                    mSum += tmpModePower;
                    if(bfTol && m%10 == 0 && tmpModePower / mSum < fTolerance) break;
                }
                nSum += mSum;
                if(bfTol && n > h && mSum / nSum < fTolerance) break;
            }
            totalSum += nSum;
            if(bfTol && nSum / totalSum < fTolerance) break;
        }

        return totalSum;
    }


    void Waveguide::OpenCSV(const std::string &filename)
    {
        const std::string d = ","; //delimeter
        bool bFileExists = bool(boost::filesystem::exists(filename));
        file.open(filename, std::ios::app);
        if(!bFileExists)
            file <<"frequency"<<d<<"magnetic_field"<<d<<"rho"<<d<<"total_power"<<d<<"TE11_power"<<"\n";

         file << std::setprecision(16);
    }

    void Waveguide::WriteCSV(const Beta &beta, const double &totalPower, const double &te11Power)
    {
        const std::string d = ","; //delimeter
        file<<beta.fc<<d;
        file<<beta.magnetic_field<<d;
        file<<beta.rho<<d;
        file<<totalPower<<d;
        file<<te11Power<<"\n";
    }

    void Waveguide::CloseCSV()
    {
        file.close();
    }

} /* namespace urchin */
