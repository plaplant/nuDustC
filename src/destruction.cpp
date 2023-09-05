/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.*/

#include "destruction.h"

#include "utilities.h"
#include "constants.h"
#include "network.h"
#include "elements.h"
#include "sputter.h"
#include "cell.h"

#include <iomanip>
#include <limits>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/format.hpp>
#include <plog/Log.h>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numeric>
#include <plog/Log.h>
#include <string>
#include <valarray>
#include <vector>
#include <sstream>

using boost::math::interpolators::makima;
using namespace boost::numeric::odeint;
using namespace std::chrono;

destruction::destruction(cell_state cell_d, network* n, params* sputARR ): 
    net (n), 
    sputARR (sputARR), 
    cell_dest (cell_d), 
    elm("data/elements.json")
{}

double destruction::Y(const double& E, const int grnid, const int gsID)
{
    using constants::pi;
    using constants::echarge_sq;
    using numbers::onehalf;
    using numbers::twothird;
    using utilities::square;
    using numbers::sp_3_441;
    using numbers::sp_2_718;
    using numbers::one;
    using numbers::sp_6_35;
    using numbers::sp_6_882;   
    using numbers::sp_1_708;
    using constants::y_pref;

    if(E<sputARR->Eth[grnid][gsID]){return 0.0;}

    // biscaro 2016 eq 7
    auto eps = sputARR->eiCoeff[grnid][gsID] * E;

    // sqrt(eps)
    auto sqrt_eps = std::sqrt(eps);

    // si(epsi) eq 6
    auto sieps = sp_3_441 * sqrt_eps * std::log(eps + sp_2_718)
                / (one + sp_6_35 * sqrt_eps + eps * (sp_6_882 * sqrt_eps - sp_1_708));
    // Eth/E
    auto eth_ratio = sputARR->Eth[grnid][gsID]/E;

    auto Si = sputARR->SiCoeff[grnid][gsID] * sieps;

    auto suffix = (one - std::pow(eth_ratio,twothird)) * square(one - eth_ratio);

    auto preret = Si * sputARR->alpha[grnid][gsID] * suffix / 
                (sputARR->u0[grnid] * (sputARR->K[grnid] * sputARR->mu[grnid][gsID] + one));

    //TODO: make sure all units are being correctly converted
    //eV_lost += 4.2E-2 * preret * sputARR->Eth[grnid][gsID];
    return y_pref * preret;
}

double destruction::Therm(const int gidx, const int gsID)
{

    double pref = sputARR->msp_2rhod[gidx]  * 
                    sqrt(sputARR->y8_piMi[gsID]*cell_dest.kT);
    auto f = [&](const double& x) { return x * exp(-x)*Y(x*cell_dest.kTeV,gidx,gsID); };
    double lowLim = sputARR->Eth[gidx][gsID]/cell_dest.kTeV;
    double Q = boost::math::quadrature::gauss_kronrod<double, 7>::integrate(f, 
                lowLim, std::numeric_limits<double>::infinity());
    return pref * Q * cell_dest.abund_moments_sizebins[gsID]; 
}

double destruction::NonTherm(const double& sidx, const int gidx, const int gsID)
{
    using constants::JtoEV;
    using constants::kB_eV;
    using constants::cm2m;
    using constants::amu2g;
    using numbers::onehalf;
    using utilities::square;
    auto iid = gidx*cell_dest.numBins+sidx;
    double x = onehalf * sputARR->miKG[gsID] * 
                square(cell_dest.vd[iid]*cm2m)*JtoEV; 
    double pref = sputARR->msp_2rhod[gidx] * cell_dest.vd[iid] * cell_dest.abund_moments_sizebins[gsID];
    return  pref * Y(x,gidx,gsID); 
}

//TODO: properly namespace `std::` (e.g. `std::pow`)
void destruction::destroy(const int gidx)
{
    using constants::pi;
    using constants::k_B;
    using constants::amu2g;
    using constants::JtoEV;
    using constants::kB_eV;
    using constants::kBeta;
    using numbers::onehalf;
    using utilities::square;
    using numbers::ten;
    for( size_t sidx =0; sidx < cell_dest.numBins; ++sidx)
    {
        double dadt = 0.0;
        auto idx = (gidx*cell_dest.numBins)+sidx;
        // remember grain sizes are in cm
        if(cell_dest.sizeBins[idx] != 0.0)
        {
            for(size_t gsID=0; gsID < cell_dest.numGas; ++gsID)
            {
                if(cell_dest.abund_moments_sizebins[gsID] == 0.0) continue;
                // s_i2 is unitless, invkT is in cgs, vd is in cm/s
                double s_i2 = sputARR->miGRAMS[gsID] * onehalf * cell_dest.invkT * square(cell_dest.vd[idx]);
                if( s_i2 > ten) 
                {
                    dadt +=  NonTherm(sidx,gidx,gsID);
                }
                else 
                {
                    dadt += Therm(gidx,gsID);
                }
            }
        }
        cell_dest.erosion_dadt[idx]+=dadt; 
        double temp_velo = calc_dvdt(cell_dest.grn_sizes[sidx],cell_dest.vd[idx],gidx) * cell_dest.dt; // in cm/s
        if(cell_dest.vd[idx] - std::abs(temp_velo) >= 0.0)
        {
            cell_dest.vd[idx] = cell_dest.vd[idx] - std::abs(temp_velo);
        }
        else{
            cell_dest.vd[idx] = 0.0;
        }
    }
}

double destruction::calc_dvdt(const double& cross_sec, const double& vd, const int grnid) // cross sec is in cm
{
    using constants::pi;
    using constants::k_B;
    using constants::amu2g;
    using numbers::eight_threeRootPi;
    using numbers::ninePi_sixyfour;
    using utilities::square;
    using numbers::one;

    double G_tot = 0;

    //TODO: dont use a reference variable for integer types
    for(size_t gsID=0; gsID < cell_dest.numGas; ++gsID)
    {
        // units of # of particles * mass in grams. might just need the mass not the * # of particles
        double m = sputARR->miGRAMS[gsID]; // should be in grams now
        double s2 = m * square(vd) / (2.*cell_dest.kT); // assumes cgs units
        G_tot += cell_dest.abund_moments_sizebins[gsID] * std::sqrt(s2) * eight_threeRootPi * 
                std::sqrt(1.+s2*ninePi_sixyfour);
    }
    auto yield = sputARR->three_2Rhod[grnid] * cell_dest.kT/(cross_sec)*G_tot;
    return  -std::abs(yield); // should be in cm/s, cross sec is in cm
}