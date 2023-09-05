/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.*/

#pragma once

#include <cmath>

namespace constants
{
const double kBoltzmann = 1.38064852E-16;
const double pressure_0 = 1.0E6;

const double pi         = M_PI;
const double shapeFactor = std::pow( 36.0 * M_PI, 1. / 3. );
const double stdP       = 1013250.0;            //# convert from dynes/cm2 to Ba, read in data should be in cgs
const double k_B        = 1.38064852E-16;        //# ergs / K
const double k_BJ       = 1.380649E-23;         //# J/K
const double istdP      = 1.0 / stdP;          //# 1/Ba
const int    N_MOMENTS  = 4;
const double ang2cm     = 1.0E-8;
const double amu2g      = 1.66054E-24;
const double amu2Kg     = 1.66054E-27;
const double JtoEV      = 6.242e+18;
const double ErgToEV    = 6.242e+11;
const double kB_eV      = 8.617333262145E-5;   // in eV
const double bohrANG    = 0.529;
const double cm2m       = 1.e-2;
const double yrs2s      = 3.1536e+7;
const double eV2K       = 8.617328149744E-5;
const double y_pref     = 4.2E-2;

// this comes from dimensional analysis
inline const double echarge_sq  = 14.399764; // [eV].[Angstrom]
inline const double echarge     = std::sqrt(echarge_sq);


inline const auto kBeta ( const double T ) 
{
    return kBoltzmann * T;
}
inline const auto equPres ( const double A, const double B, const double T ) 
{
    return std::exp ( -A / T + B );
}

} // namespace constants

namespace numbers
{
    inline const double one     = 1.0;
    inline const double two     = 2.0;
    inline const double three   = 3.0;
    inline const double four    = 4.0;
    inline const double ten     = 10.0;
    inline const double eight   = 8.0;
    inline const double nine    = 9.0;
    inline const double sixtyfour = 64.0;

    inline const double onehalf     = one / two;
    inline const double onethird    = one / three;
    inline const double onefourth   = one / four;
    inline const double twothird    = two / three;

    // constants for calculating biscaro 2016 
    inline const double eight_threeRootPi   = eight / (three * std::pow(M_PI,onehalf));
    inline const double ninePi_sixyfour     = nine * M_PI / sixtyfour;
    inline const double sp_3_441   = 3.441;
    inline const double sp_2_718   = 2.718;
    inline const double sp_6_35    = 6.35;
    inline const double sp_6_882   = 6.882;   
    inline const double sp_1_708   = 1.708;
} // namespace numbers