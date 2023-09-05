/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.*/

#include "nucleation.h"

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

nucleation::nucleation ( cell_state cell_n, network* n ): net (n), cell_nuc (cell_n), elm("data/elements.json")
{}

void
nucleation::nucleate(const std::vector<double>& x, const double time)
{
  using constants::kBeta;
  using constants::equPres;
  using constants::pi;
  using constants::istdP;
  using constants::k_B;
  using constants::stdP;
  using constants::amu2g;
  using constants::kB_eV;

  std::fill(std ::begin(cell_nuc.erosion_dadt), std ::end(cell_nuc.erosion_dadt), 0.0);

  for (size_t i = 0; i < cell_nuc.numReact; ++i) {
    auto reaction_idx = net->nucleation_reactions_idx[i];
    auto num_ks       = net->ks_lists_idx[i].size();
    auto key_spec_idx = net->ks_lists_idx[i][0];
    if (num_ks > 1) {
      for (size_t kidx = 0; kidx < num_ks; ++kidx) {
        auto r_idx = net->ks_lists_idx[i][kidx];
        if (x[key_spec_idx] > x[r_idx])
          key_spec_idx = r_idx;
      }
    }
    cell_nuc.parts[i].ks_idx = key_spec_idx;
    if (x[key_spec_idx] < CELL_MINIMUM_ABUNDANCE) {
      cell_nuc.parts[i].is_nucleating = 0;
      cell_nuc.S[i]                   = 0.0;
      cell_nuc.Js[i]                  = 0.0;
      cell_nuc.dadt[i]                = 0.0;
      cell_nuc.ncrit[i]               = 0.0;
      continue;
    }

    // std::string elem_sym = net->species[cell_nuc.parts[i].ks_idx];
    // auto at_mass = elm.elements.at(elem_sym).mass;
    // cell_nuc.parts[i].ks_react_mass = at_mass * amu2g;
    cell_nuc.parts[i].ks_react_mass =
      elm.elements.at(net->species[cell_nuc.parts[i].ks_idx]).mass * amu2g;

    std::vector<double> react_nu;
    std::vector<int> react_idx;
    auto stoich_ks = 0.0;
    for (const auto& kv: net->nucleation_species_count[i]) {
      if (kv.first == key_spec_idx) {
        stoich_ks = kv.second;
      }
      react_idx.emplace_back(kv.first);
    }
    for (const auto& kv: net->nucleation_species_count[i])
      react_nu.emplace_back(kv.second / stoich_ks);
    cell_nuc.parts[i].react_idx = react_idx;
    cell_nuc.parts[i].react_nu  = react_nu;

    double psum = 0.0;
    for (size_t ridx = 0; ridx < react_idx.size(); ridx++) {
      if (x[react_idx[ridx]] != 0) {
        if (react_idx[ridx] != cell_nuc.parts[i].ks_idx)
          psum = psum + log(x[react_idx[ridx]] * cell_nuc.kT * istdP) * react_nu[ridx];
      }
    }

    double c1 = x[key_spec_idx];
    cell_nuc.cbars[i] = cell_nuc.init_abund[key_spec_idx] * cell_nuc.volume_0 / cell_nuc.volume;

    auto delg_reduced = (net->reactions[reaction_idx].alpha / cell_nuc.temperature -
                         net->reactions[reaction_idx].beta) + psum;
    cell_nuc.parts[i].lnS  = log(c1 * cell_nuc.kT * istdP) + delg_reduced;
    cell_nuc.parts[i].catS = (stdP / cell_nuc.kT) * exp(-delg_reduced);

    double w = 1.0;
    for (size_t ridx = 0; ridx < react_idx.size(); ++ridx) 
    {
      if (react_idx[ridx] != cell_nuc.parts[i].ks_idx)
        w = w + react_nu[ridx];
    }
    double Pii = 1.0;
    for (size_t ridx = 0; ridx < react_idx.size(); ++ridx) 
    {
      if (react_idx[ridx] != cell_nuc.parts[i].ks_idx)
        Pii = Pii * pow(x[react_idx[ridx]] / c1, react_nu[ridx]);
    }
    if (cell_nuc.parts[i].lnS > 0.0) {
      double iw = 1. / w;
      Pii       = pow(Pii, iw);
      double mu = 4.0 * pi * pow(net->reactions[reaction_idx].a_rad, 2.) *
                  net->reactions[reaction_idx].sigma / cell_nuc.kT;
      double expJ = -4.0 / 27.0 * pow(mu, 3.) / pow(cell_nuc.parts[i].lnS, 2.);
      double Jkin = pow(2.0 * net->reactions[reaction_idx].sigma /
                          (pi * cell_nuc.parts[i].ks_react_mass),0.5);
      cell_nuc.S[i]    = exp(cell_nuc.parts[i].lnS);
      cell_nuc.Js[i]   = net->reactions[reaction_idx].omega0 * Jkin * c1 * c1 * Pii * exp(expJ);
      cell_nuc.dadt[i] = net->reactions[reaction_idx].omega0 *
                     pow(0.5 * cell_nuc.kT / (pi * cell_nuc.parts[i].ks_react_mass), 0.5) *
                     c1 * (1. - 1. / cell_nuc.S[i]);
      cell_nuc.ncrit[i] = pow(2.0 / 3.0 * (mu / cell_nuc.parts[i].lnS), 3.0) + iw;


      auto momIDX = cell_nuc.numSpec + constants::N_MOMENTS * i;
    } 
    else 
    { // we want to force these to zero in case a value is unchanged for
      // the next timestep
      cell_nuc.S[i]     = 0.0;
      cell_nuc.Js[i]    = 0.0;
      cell_nuc.dadt[i]  = 0.0;
      cell_nuc.ncrit[i] = 0.0;
    }
    if (cell_nuc.ncrit[i] > 0.0) 
    {
      cell_nuc.parts[i].grains_nucleating = cell_nuc.Js[i] * cell_nuc.ncrit[i];
      cell_nuc.parts[i].is_nucleating     = 1;
    } 
    else 
    {
      cell_nuc.parts[i].grains_nucleating = 0.0;
      cell_nuc.parts[i].is_nucleating     = 0;
    }
  }
}