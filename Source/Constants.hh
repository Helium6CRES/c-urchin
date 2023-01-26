#ifndef CONSTANTS_HH_
#define CONSTANTS_HH_

/**
 * This class contains various fundamental constants.
 * Values are taken from PDG edition 2006, unless pointed out otherwise. The naming conventions are: normal name for SI units, a suffix _unit for something else.
 **/

class Constants
{
public:
    Constants() = delete;

    constexpr static double Pi()
    {
        return 3.141592653589793238462643383279502884L;
    } //!< pi

    constexpr static double C()
    {
        return 299792458.0;
    } //!< c im m/s

   constexpr static double Q()
    {
        return 1.60217653E-19;
    } //!< elementary charge  in C(>0)

    //EM coupling constants
    constexpr static double EpsNull()
    {
        return 8.854187817E-12;
    } //!< epsilon0, Constant of Newtons force.

    constexpr static double MuNull()
    {
        return 4.E-7 * Pi();
    }//!< permeability of free space

    //masses
    constexpr static double M_el_kg()
    {
        return 9.1093826E-31;
    } //!< electron mass in kg

    constexpr static double M_el_eV()
    {
        return 510.998918E3;
    } //!< electron mass in ev/c^2
};
#endif 
