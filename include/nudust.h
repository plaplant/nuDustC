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
#include "network.h"
#include "cell.h"
#include "sputter.h"
#include "sput_params.h"

#include <vector>
#include <map>
#include <string>

namespace options = boost::program_options;

class nuDust
{
  sput::sputter_list_t            sputter;
  std::vector<std::string>        initial_elements;
  std::map<uint32_t, cell_input>  cell_inputs;
  std::vector<cell>               cells;
  std::vector<cell>             RScells;
  std::map<uint32_t, cell_input> RScell_input;
  network               net;
  params                sputARR;
  std::vector<double>   init_size_bins; // sizes in the size dist
  std::vector<double>   init_bin_edges; // the edges of the size dist
  std::vector<double>   emptySD_VD;     // empty array of length numBins*num_grains for the SD and ve
  std::vector<double>   size_bins_init; // user defined initial size distribution
  
  std::string name;
  std::string nameRS;

  int par_size, par_rank;
  int numBins;

public:
  configuration nu_config;
  nuDust(const std::string& config_filename, int sz, int rk);
  nuDust(const std::string& config_filename);
  virtual ~nuDust() {}

  void load_network();
  void load_initial_abundances();
  void load_sizeDist();
  void load_shock_params();
  void load_sputter_params();
  void find_shock();
  void load_environment_data();
  void load_outputFL_names();
  void gen_size_dist();
  void account_for_pileUp();
  void gen_shock_array_frm_val();
  void create_restart_cells(int cid);
  void generate_sol_vector();
  void create_simulation_cells();
  int get_element_index(const std::string& elem) const;
  void premake(const int s1, const int s2, const int sp, const int cell_id);
  void run();
  void cleanup() {}
};
