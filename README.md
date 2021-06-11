# Self Avoiding Random Walk

## Overview

This Self-Avoiding Random Walk (SARW) generates non-overlapping and randomly arranged polymer chains for molecular dynamic simulations.

Instead of creating fully atomistic polymer chain, this method uses united atoms to model each monomer, where all hydrogens are hidden. For example, CH3 and CH2 molecules are represented by particles of mass 15 and 14 respectively. The polymer chain

This method creates a chain of particles at a fixed distance equal to the carbon-carbon bond length of the polymer being represented, and at a fixed carbon-carbon-carbon bond angle. The orientation of the bond angle is randomly determined.

The method is implemented in C++ and is contained in the following files:
- `sarw.cpp`: the main file with all key subroutines to create the polymer chains
- `sarw.h`: the header file for the above with definitions for key variables and data structures
- `InputMap.cpp`: program to handle the input file with user-defined variables
- `InputMap.h`: header file for the above
- `PeriodicLookup.h`: creates a lookup table for finding potential overlapping particles. This is responsible for significant speedup.

## How to run

- list all the necessary user-defined variables and their desired values in an ASCII. Eg, `pe.in`
- build the `sarw` binary using CMake
- Use the binary and input file to create polymer chains. Eg, `sarw -i pe.in`

## Description of user-defined variables

<!-- All model parameters are defined in the sarw.h file -->

|Parameter | Description | Possible values | Default |
|---|---|---|---|
|`boundary`|Turn off periodic boundary conditions on one of the cartesian directions (similar to LAMMPS syntax).|3 options: ppf, pfp, fpp|ppp|
|`boxSize`|Size of the cuboid shaped simulation box in Angstroms|Any set of 3 integers separated by commas|10,10,10|
|`graftFraction`|Fraction of total chains which need to be grafted at locations listed in 'grafts.dat'|Any fraction less than 1|0.0|
|`growthBias`|A vector used to grow the chains in specific cartesian direction.|Any set of 3 integers separated by commas. Each bias can either be 0 +1 or -1|0,0,0|
|`logProgress`|Whether to show overall progress of chain growth in the output file.|yes/no|yes|
|`logSteps`|Whether to show details about every Monte Carlo attempt, and each function call to grow the chains.|yes/no|no|
|`maxTrials`|Set an upper limit on the Monte Carlo trials.|Any integer|100|
|`minDist`|Set a minimum distance allowed between non-bonded united-atoms.|Any positive floating number|2.0|
|`nChains`|Set the number of chains that are grown in the simulation box. Some of them maybe grafted if `graftFraction` is non-zero.|Any integer|1|
|`polymer`|Name of the polymer used in the header of LAMMPS data file|Any string|polymer|
|`roundRobin`|Whether the chains should be grown in round-robin fashion or randomly|yes/no|yes|
|`targetMassDensity`|Set the mass density of simulation box in gram per cubic centimeter|Any positive fraction|0.92|

## Outputs
- `polymer.XYZ`: type and cartesian coordinates of all particles, can be opened in Ovito and XMD
- `polymer.XSF`: box size, and coordinates of all particles, particles are not distinguishable as Carbon is used to represent both the types, can be opened in VESTA
- `polymer.data`: data file in the prescribed LAMMPS format, can be read in using 'read_data' command in LAMMPS,

## TODO

- [ ] arbitrary constraints on a particular component of positions (eg. grow polymer backbone on a plane)
- [x] Use Classes/Objects to pass the data around
- [x] Input file for reading the customizable parameters
- [x] Periodic Lookup table
- [ ] Construct backbone and then add sideChains
- [ ] Total time taken
