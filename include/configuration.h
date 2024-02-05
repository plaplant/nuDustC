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

#include <boost/program_options.hpp>
#include <string>

namespace options = boost::program_options;

struct configuration
{
  double ode_dt_0;
  double ode_dt_min;
  double ode_dt_max;
  double ode_abs_err;
  double ode_rel_err;
  double low_sd_exp;
  double high_sd_exp;
  double shock_velo;
  double shock_temp;
  double pile_up_factor;
  double sim_start_time;


  int io_dump_n_steps;
  int io_restart_n_steps;
  int bin_number;

  int do_destruction;
  int do_nucleation;

  std::string network_file;
  std::string sizeDist_file;
  std::string abundance_file;
  std::string shock_file;
  std::string environment_file;

  // used to differentiate runs or models
  std::string mod_number;
  

  options::options_description desc;
  options::variables_map vm;

  configuration();
  void read_config(const std::string& filename);
};
