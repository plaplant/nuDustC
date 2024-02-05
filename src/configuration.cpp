/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.*/

#include "configuration.h"

#include <boost/program_options.hpp>
#include <fstream>
#include <plog/Log.h>
#include <string>
#include <math.h>

configuration::configuration() : desc ( "configuration" )
{

    // ODE description
    desc.add_options() ( "ode_dt_0", options::value<double> ( &ode_dt_0 )->default_value ( 1.0E-3 ), "initial time step" );
    desc.add_options() ( "ode_abs_err", options::value<double> ( &ode_abs_err )->default_value ( 1.0E-6 ), "solver absolute error criteria" );
    desc.add_options() ( "ode_rel_err", options::value<double> ( &ode_rel_err )->default_value ( 1.0E-6 ), "solver relative error criteria" );
    desc.add_options() ( "ode_dt_min", options::value<double> ( &ode_dt_min )->default_value ( 1.0E-6 ), "solver minimum allowed dt" );
    desc.add_options() ( "ode_dt_max", options::value<double> ( &ode_dt_max )->default_value ( 1.0E2 ), "solver max allowed dt" );
    
    // Input data files
    desc.add_options() ( "sizeDist_file", options::value<std::string> ( &sizeDist_file ), "file with size distributions" );
    desc.add_options() ( "environment_file",options::value<std::string>(&environment_file),"file with stellar environment variables");
    desc.add_options() ( "network_file", options::value<std::string> ( &network_file ), "file with network" );
    desc.add_options() ( "abundance_file", options::value<std::string> ( &abundance_file ), "file with inital abundances" );
    desc.add_options() ( "shock_file", options::value<std::string> ( &shock_file ), "file with inital shock parameters per cell" );

    // describe dust size distribution for binning
    desc.add_options() ( "size_dist_min_rad_exponent_cm", options::value<double> ( &low_sd_exp )->default_value ( NAN ), "exponent for the min radius of the size dist. in cm" );
    desc.add_options() ( "size_dist_max_rad_exponent_cm", options::value<double> ( &high_sd_exp )->default_value ( NAN ), "exponent for the max radius of the size dist. in cm" );
    desc.add_options() ( "number_of_size_bins", options::value<int> ( &bin_number )->default_value ( NAN ), "bin number" );
    
    // determine if doing nucleation and/or destruction
    desc.add_options() ( "do_destruction", options::value<int> ( &do_destruction )->default_value (0), "do destruction calculations" );
    desc.add_options() ( "do_nucleation", options::value<int> ( &do_nucleation )->default_value (0), "do nucleation calculations" );

    // print out and save to file controls
    desc.add_options() ( "io_restart_n_steps",options::value<int>(&io_restart_n_steps)->default_value(1000),"write restart file to disk every n steps");
    desc.add_options() ( "io_dump_n_steps",options::value<int>(&io_dump_n_steps)->default_value(1000), "write dump file to disk  every n steps");
    

    // user specified shock parameters. shock parameters are applied to all cells regardless of depth
    desc.add_options() ( "pile_up_factor", options::value<double> ( &pile_up_factor )->default_value ( 1.0 ), "shock pile up factor" );
    desc.add_options() ( "shock_velo", options::value<double> ( &shock_velo )->default_value ( NAN ), "shock velocity" );
    desc.add_options() ( "shock_temp", options::value<double> ( &shock_temp )->default_value ( NAN ), "shock temperature" );

    // specify model start time. used to adjust the abundances for any given time. assumes homologous expansion
    desc.add_options() ( "sim_start_time", options::value<double> ( &sim_start_time )->default_value ( NAN ), "simulation start time" );

    // remnant for models. can set to anything
    desc.add_options() ( "mod_number", options::value<std::string> ( &mod_number ), "model number" );
}

// reading in the config file and declaring variables
void
configuration::read_config(const std::string& config_filename)
{
  std::ifstream config_file(config_filename.c_str());
  vm = options::variables_map();
  options::store(options::parse_config_file(config_file, desc), vm);
  options::notify(vm);
}
