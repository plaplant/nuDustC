

# nuDustC
Nucleating Dust Code in C++

***All Units are in CGS***

### To Build & Run
**Required:** OpenMP, MPI, Boost, PLOG.
This build uses cmake. 

### Plog 
Plog is a header only package. Running cmake fetches Plog, but it can also be aquired by:

'''
git clone git@github.com:SergiusTheBest/plog.git
'''


# Building nuDustC++
In order to build nuDustC++, go to the head of the git repository (nudustc/) and run:

```
$> mkdir build;
$> cd build;
$> cmake .. ;
$> make;
```

Run nuDustc with

```
$> ./nudust++ -c data/inputs/test_config.ini
```

You may need to add the paths of the packages to '''LD_LIBRARY_PATH'''. If using conda:

```
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:~/miniconda3/lib
```

nuDustC++ can then be run by:

```
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:~/miniconda3/lib ./nudust++ -c data/inputs/test_config.ini
```

### Inputs
Required: Config file. This lists the various input information such as data files, integration parameters, and calculation options.

# Configuration file

### Integration Parameters
*ode_dt_0*: Initial timestep for integration.

*ode_abs_err*: The Solver's absolute error.

*ode_rel_err*: The Solver's relative error.

*ode_dt_min*: The minimum allowed timestep.

*de_dt_max*: The maximum allowed timestep.

    
### Data Files
*sizeDist_file*: This describes the size distribution for the model. Each cell has one line. Each line is an array of size distribtuions of grain species in the order specified in the header line.

*environment_file*: This contains the trajectory data for each timestep. The time is specified on a single line. Below, each cell's is describe in a single line: cell_ID, temperature (K), volume(cm^3), density(g/cm^3), pressure (Ba), velocity (cm/s), radius (cm).

*network_file*: This includes the chemical network of grain reactions. Each grain species takes up one line in order: reactants, "->", products, "|", key species, Gibb's free energy 'A' term (A/10^4 K), Gibb's free energy 'B' term, surface energy of the condensate (ergs/cm^2), radius of condensate (angstroms). 

*abundance_file*: This lists the names of gas species in the header. Each cell has one line listing: cell ID, number density for each gas species. 

*shock_file*: This contains information on a shock. Each cell has one line: cell ID, the time of the shock, the shock temperature, the shock velocity.

### Size Distribution Parameters
*size_dist_min_rad_exponent_cm*: The exponent of the left edge of the distribtuion.

*size_dist_max_rad_exponent_cm*: The exponent of the right edge of the distribtion.

*number_of_size_bins*: The number of size bins.
    
### Control Nucleation and Destruction
*do_destruction*: 1, if doing destruction, 0 if no destruction.

*do_nucleation*: 1, if soind nucleation, 0 if no nucleation.

If both are set to '1', both destruction and nucleation are calculated. 

### Data Output Controls
*io_disk_n_steps* : Number of cycles until a dump file is updated.  

*io_screen_n_steps*: Number of cycles until a restart file is updated. 

### User Specified Shock Parameters
*pile_up_factor*: This is used to calculate the increase in density when a shock passes through. The density if multiplied by this number. 

*shock_velo*: This is the velocity of the shock that will be applied to all cells. 

*shock_temp*: This is the temperature of the shock that will be applied to all cells.

*sim_start_time*: This is the time of the shock. It is used to calculate the change in volume, assuming homologous expansion.

### Nomenclature
*mod_number*: This specifies the model used in order to differentiate multiple runs.

# Functionality
Nucleation & Destruction
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density, Pressure, Velocity), Abundance File, & Network File

Nucleation
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density), Abundance File, & Network File

Destruction With Shock Times and Velocities from a file
  Required Input Files: Shock File (Cell #, Time, Shock Temperature, Shock Velocity), Pile up factor, Size Distribution File, Abundance File, & Network File

Destruction With User Input Shock Temperature & Velocity
  Required Input Files: Shock Velocity, Shock Temperature, Shock Time, Pile up factor, Size Distribution File, Abundance File, & Network File
