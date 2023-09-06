
# nuDustC
Nucleating Dust Code in C++

***All Units are in CGS***

### To Build & Run
**Required:** OpenMP, MPI, Boost
This build uses cmake. Run
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
### Inputs
Required: Config file. This lists the various input information such as data files, integration parameters, and calculation options.

### Functionality
Nucleation & Destruction
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density, Pressure, Velocity), Abundance File, & Network File

Nucleation
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density), Abundance File, & Network File

Destruction With Shock Times and Velocities from a file
  Required Input Files: Shock File (Cell #, Time, Shock Temperature, Shock Velocity), Pile up factor, Size Distribution File, Abundance File, & Network File

Destruction With User Input Shock Temperature & Velocity
  Required Input Files: Shock Velocity, Shock Temperature, Shock Time, Pile up factor, Size Distribution File, Abundance File, & Network File
