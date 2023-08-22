# nuDustC
Nucleating Dust Code in C++

**All Units are in CGS**


# Inputs
Required: Config file. This lists the various input information such as data files, integration parameters, and calculation options.

# Functionality
Nucleation & Destruction In Time
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density, Pressure, Velocity), Abundance File, & Network File

Nucleation In Time
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density), Abundance File, & Network File

Destruction In Time
  Required Input Files: Hydrodynamical Trajectory file (Time, Temperature, Volumes, Density), Size Distribution File, Abundance File, & Network File

Destruction With Shock Times Specified
  Required Input Files: Shock File (Cell #, Time, Shock Temperature, Shock Velocity), Size Distribution File, Abundance File, & Network File

Destruction With User Input Shock Temperature & Velocity
  Required Input Files: Shock Velocity, Shock Temperature, Pile up factor, Size Distribution File, Abundance File, & Network File
