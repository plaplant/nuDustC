

# nuDustC
Nucleating Dust Code in C++

***All Units are in CGS***

### To Build & Run
**Required:** OpenMP, MPI, Boost, SunDials, Plog.
This build uses cmake. 

### Plog 
Plog is a header only package. Running cmake fetches Plog, but it can also be aquired by:

```
git clone git@github.com:SergiusTheBest/plog.git
```

### Boost
If using conda, install boost as:

```
conda install conda-forge::boost
```

Further documentation and installation of Boost can be found at:

https://www.boost.org/

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

You may need to add the paths of the packages to *LD_LIBRARY_PATH*. If using conda nuDustC++ can then be run by:

```
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:~/miniconda3/lib ./nudust++ -c data/inputs/test_config.ini
```

Or specifying the paths to the packages and libraries.

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
There are 5 main calculation paths including destruction or nucleation calculations. Destruction has 2 main branches based on if the user is specifying shock values or reading them in from a file. 

Nucleation & Destruction With User Input Shock Temperature & Velocity
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density, Pressure, Velocity), Abundance File, & Network File, Shock Velocity, Shock Temperature, Shock Time, Pile up factor, Size Distribution File or Size Parameters

Nucleation & Destruction With Shock Times and Velocities from a file
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density, Pressure, Velocity), Abundance File, & Network File, Shock Velocity, Shock Temperature, Shock Time, Pile up factor, Size Distribution File or Size Parameters

Nucleation
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density), Abundance File, & Network File, Size Distribution File or Size Parameters

Destruction With Shock Times and Velocities from a file
  Required Input Files: Shock File (Cell #, Time, Shock Temperature, Shock Velocity), Size Distribution File or Size Parameters, Abundance File, & Network File

Destruction With User Input Shock Temperature & Velocity
  Required Input Files: Shock Velocity, Shock Temperature, Shock Time, Pile up factor, Size Distribution File or Size Parameters, Abundance File, & Network File

# Selecting Integrators and Interpolators
The integrator is setup in *src/cell.cpp* in the *solve()* funtion. nuDustC++ comes defaulted with a runge-kutta doPri 5 integrator. Additional information on the available integrators offered by Boost can be found at:

https://www.boost.org/doc/libs/1_78_0/libs/numeric/odeint/doc/html/index.html

nuDustC++ currently uses a makima 1-D interpolator. The interpolator is defined in *include/cell.h* in the 'cell' class declaration. Additional interpolators offered by Boost can be found at:

https://www.boost.org/doc/libs/1_78_0/libs/math/doc/html/interpolation.html

# Testing
To run a test of nudustc++,

```
./nudustc++ -c data/inputs/test_config.ini
```

This will produce data files in the build directory's "output/" directory and restart data in the "restart/" directory.

# Restarting a Run
nuDustC++ automatically checks for restart files when creating each cell. If a restart file is found, that data is loaded into the cell object. If no restart file is found, the cell is initialized with data from the input files. Make sure the same config file used to start the run is selected when restarting. 

# Common Pitfalls
If the compiler cannot find required packages or libraries, make sure LD_LIBRARY_PATH is up to date and points to the location of each package or library.

Ensure the configuration file points to the location of each input file and contains the needed parameters for the calculation path. 


