#include <iostream>
#include <math.h>
#include <string>
#include <boost/program_options.hpp>

#include "Beta.hh"
#include "CircularWaveguide.hh"
#include "RectangularWaveguide.hh"

namespace po = boost::program_options;

int main(int ac, char* av[])
{
    unsigned N, M, H;
    unsigned nFrequencies, nRhos, nFields;
    double startFrequency, startRho, startField;
    double endFrequency, endRho, endField;
    double fTol;
    bool bCircularWaveguide, bTE, bTM;
    std::string outputFilename;

    po::options_description desc("Allowed options");
    try {
        desc.add_options()
            ("help,-h", "Produce help message")
            ("max-N-mode,N", po::value<unsigned>(&N)->default_value(1600), "Maximum possible mode value in summing over n")
            ("max-M-mode,M", po::value<unsigned>(&M)->default_value(1600), "Maximum possible mode value in summing over m")
            ("max-harmonic,H", po::value<unsigned>(&H)->default_value(200), "Maximum possible harmonic h")
            ("tolerance,t", po::value<double>(&fTol)->default_value(0), "Truncate h-sum if sum changes <fTol (0 for no skipping)")
            ("output-filename,o", po::value<std::string>(&outputFilename)->default_value("out.csv"), "Filename for simulation outputs")
            ("circular-waveguide,c", po::value<bool>(&bCircularWaveguide)->default_value(true), "Circular Waveguide? False for rectangular.")
            ("TE,e", po::value<bool>(&bTE)->default_value(true), "Add TE Modes to sum")
            ("TM,m", po::value<bool>(&bTM)->default_value(false), "Add TM Modes to sum")
            ("number-frequencies", po::value<unsigned>(&nFrequencies)->default_value(1), "Number of frequencies")
            ("number-rhos", po::value<unsigned>(&nRhos)->default_value(1), "Number of rho values")
            ("number-fields", po::value<unsigned>(&nFields)->default_value(1), "Number of field values")
            ("start-frequency", po::value<double>(&startFrequency)->default_value(19e9), "Minimum frequency")
            ("start-rho", po::value<double>(&startRho)->default_value(0.), "Minimum rho value")
            ("start-field", po::value<double>(&startField)->default_value(2.), "Minimum B field value")
            ("end-frequency", po::value<double>(&endFrequency)->default_value(19e9), "Maximum frequency")
            ("end-rho", po::value<double>(&endRho)->default_value(0.), "Maximum rho value")
            ("end-field", po::value<double>(&endField)->default_value(2.), "Maximum B field value")
        ;

        po::variables_map vm;        
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 0;
        }
    }
    catch(std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        std::cout << desc << "\n";
        return 1;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }
    if(!bTE && !bTM)
    {
        std::cerr << "Neither TE nor TM set: Returning!"<<std::endl;
        return 2;
    }

    auto dX = [](double a, double b, unsigned N ){ return N==1 ? 0 : (b-a)/(N-1);};

    const double dB = dX(startField, endField, nFields);
    const double dFrequency = dX(startFrequency, endFrequency, nFrequencies);
    const double dRho = dX(startRho, endRho, nRhos);

    const double aRadius = 0.00578;
    const unsigned nRoots = 1600;

    urchin::Waveguide *w;
    if(bCircularWaveguide)
    {
        w = new urchin::CircularWaveguide(aRadius,nRoots);
    }
    else
    {
        w = new urchin::RectangularWaveguide(10.668e-3, 4.318e-3); //Hardcoded WR42 waveguide
    }

    w->OpenCSV(outputFilename);

    //Loop over betas to simulate
    double Ptot;
    double B = startField;
    
    for(unsigned nb = 0; nb < nFields; ++nb)
    {
        double fc = startFrequency;
        for(unsigned nf = 0; nf < nFrequencies; ++nf)
        {
            double rho = startRho;
            for(unsigned nr = 0; nr < nRhos; ++nr)
            {
                urchin::Beta b(fc,B,rho);
                Ptot = 0;
                if(bTE) Ptot += w->TotalPower(N,M,H,fTol,b,true);
                if(bTM) Ptot += w->TotalPower(N,M,H,fTol,b,false);
                double PTE11 = w->ModePower(1,1,1,b,true);
                w->WriteCSV(b,Ptot,PTE11);

                rho += dRho;
            }
            fc +=  dFrequency;
        }
        B += dB;
    }

    w->CloseCSV(); //perhaps should be in destructor?
    delete w;

    std::cout<<"Done!"<<std::endl;

    return 0;
}



    
