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

#include "configuration.h"
#include "elements.h"
#include "network.h"
#include "sput_params.h"
#include "sputter.h"
#include "utilities.h"
#include "constants.h"

#include <vector>
#include <string>
#include <map>
#include <chrono>
#include <fstream>
#include <boost/format.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/serialization/array_wrapper.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

using boost::math::interpolators::makima;

const size_t CELL_MAX_STEPS              = 10000000000000;
const size_t CELL_MAXIMUM_STEPPER_RESETS = 1000;
const double CELL_MINIMUM_ABUNDANCE      = 1.0E-1;

const double NO_REBINNING_MAX = 1e4;
typedef std::vector<double> abundance_v;
const double min_dadt = 1e-12;

struct cell_input
{
  int                 inp_norebin_ntimes = 0;

  std::vector<double> inp_abund;
  std::vector<double> inp_init_abund;

  // from enviornment file
  std::vector<double> inp_times;
  std::vector<double> inp_temp;
  std::vector<double> inp_volumes;
  std::vector<double> inp_rho;
  std::vector<double> inp_pressure;
  std::vector<double> inp_velo;
  std::vector<double> inp_x_cm;

  // nucleation
  std::vector<double> inp_dadt;

  // destruction 
  std::vector<double> inp_erosion_dadt;
  std::vector<double> inp_destBins;
  std::vector<double> inp_vd;

  // this is from the user specified config params or from the shock params file
  double inp_shock_time;
  double inp_shock_temp;

  // size distribution
  // can be calculated from the user defined config or taken from the size distribution file
  std::vector<double> inp_binSizes;
  std::vector<double> inp_binEdges;
  std::vector<double> inp_size_dist;

  // user defined start time. assumes homologous expansion to compute from cell time. 
  // only use for destruction, no nucleation
  double sim_start_time; // this is from the config file
  double inp_cell_time; // this is read in from the size dist file

  // for destruction calculate shock values from enviornment file
  std::vector<double> inp_shock_times_arr;
  std::vector<double> inp_shock_velo_arr;
  std::vector<bool> inp_shock_bool_arr;

};

struct cell_partial
{
  double reaction_idx;
  double saturation = 0.0;

  bool is_nucleating     = 1;
  double critical_size   = 0.0;
  double nucleation_rate = 0.0;

  double grains_nucleating;

  // sms added this here, may move
  double dadt = 0;
  int ks_idx;
  double ks_react_mass;
  double cbar;
  std::vector<double> r_nu;
  double catS = 0.0;
  double lnS  = 0.0;
  std::vector<double> react_nu;
  std::vector<int> react_idx;
};

struct cell_state
{
  double temperature;
  double volume;
  double volume_0;
  double dvolume;
  double rho;
  double drho;
  double pressure;
  double dP;

  double kT;
  double kTeV;
  double invkT;
  double eV_lost;

  std::vector<double> init_abund;
  std::vector<double> cbars;
  std::vector<double> dadt;
  std::vector<double> Js;
  std::vector<double> ncrit;
  std::vector<double> S;
  //std::vector<double> g_change;
  std::vector<double> abund_moments_sizebins;

  std::vector<cell_partial> parts;

  std::vector<double> grn_sizes;
  std::vector<double> edges;
  double start_time;
  double time;
  double dt;

  std::vector<double> vd;
  std::vector<double> sizeBins;
  std::vector<double> runningTot_size_change;

  int numReact;
  int numBins;
  int sd_len;
  int numGas;
  int numSpec; 

  int len_abund_and_mom;

  std::vector<double> rebin_chng;
};

class cell
{
public:
  network*        net;
  params*         sputARR;
  configuration*  config;
  uint32_t        cid;
  xkin::element_list_t elm;
  std::vector<cell_state> solution_states;
  std::vector<double> init_SD;
  std::vector<bool> reaction_switch;

  std::vector<double> env_times;
  std::vector<double> env_temp;
  std::vector<double> env_volumes;
  std::vector<double> env_rho;
  std::vector<double> env_pressure;

  std::vector<double> env_shock_times;
  std::vector<double> env_shock_velo;
  std::vector<double> env_shock_bool;

  // define interpolator. probably not the best method but it works.
  std::vector<double> fake1{ 1, 2, 3, 4 }, fake2{ 1, 4, 8, 16 };
  boost::math::interpolators::makima<std::vector<double>> env_temp_interp = makima(std::move(fake1), std::move(fake2));
  std::vector<double> fake3{ 1, 5, 9, 12 }, fake4{ 1, 4, 8, 16 };
  boost::math::interpolators::makima<std::vector<double>> env_volume_interp =
    makima(std::move(fake3), std::move(fake4));
  std::vector<double> fake7{ 1, 2, 3, 4 }, fake8{ 1, 4, 8, 16 };
  boost::math::interpolators::makima<std::vector<double>> env_rho_interp =
    makima(std::move(fake7), std::move(fake8));
  std::vector<double> fake11{ 1, 2, 3, 4 }, fake12{ 1, 4, 8, 16 };
  boost::math::interpolators::makima<std::vector<double>> env_pressure_interp =
    makima(std::move(fake11), std::move(fake12));
  std::vector<double> fake13{ 1, 2, 3, 4 }, fake14{ 0, 1, 1, 0 };
  boost::math::interpolators::makima<std::vector<double>> env_shock_bool_interp =
    makima(std::move(fake13), std::move(fake14));
  std::vector<double> fake15{ 1, 2, 3, 4 }, fake16{ 0, 1, 1, 0 };
  boost::math::interpolators::makima<std::vector<double>> env_shock_velo_interp =
    makima(std::move(fake15), std::move(fake16));

  cell_state cell_st;

  bool integration_abandoned;

  void check_reactions(const std::vector<double>& x);
  bool check_solution(const std::vector<double>& x);
  void rebin (const std::vector<double>& x, std::vector<double>& dxdt);
  void calc_state_vars(const std::vector<double>& x, const double time);
  void nucleate(const std::vector<double>& x);
  void destroy();
  void add_new_grn(const std::vector<double>& x);
  double calc_dvdt(const double& cross_sec, const double& vd, const int grnid);
  double Y(const double& E, const int grnid, const int gasid);
  double Therm(const int grnid, const int gasid);
  double NonTherm(const double& sidx, const int grnid, const int gsID);
  void set_init_data(const spec_v& init_s, const cell_input& init_data);
  void set_env_data(const cell_input& input_data);

public:
  cell(network* n,
       params* sputArr,
       configuration* con,
       uint32_t id,
       const spec_v& init_s,
       const cell_input& input_data);

  virtual ~cell(){};
  void solve();
  void operator()(const std::vector<double>& x, std::vector<double>& dxdt, const double t);

};
