# Overview

Using Self-Avoiding Random Walk (SARW) to generate polymer chains for MD simulations.
Each monomer is modelled by a united atom, where all hydrogens are contracted into the carbon atoms.  For example, CH3 and CH2 are modelled as single particle of different kinds, CH3 is a particle of mass 15 a.u. and CH2 has mass of 14 a.u.
Each subsequent united atom is placed randomly on a cone of half-angle = (180 - tetrahedral angle)=75.5 degress to capture the sp3 bond angle.

## Description of Model Parameters

All model parameters are in the main.h file

- nChains = number of chains. These many number of seeds are randomly placed in the allowable region of simulation boxx

- minChainLength = target length of each chain. Propagation of chains do not stop when they reach this length.

- mass_density = desired mass density of polymer (g/cc)

- minDist = minimum distance defined between non-bonded monomers

- maxTrials = upper limit on the MC steps

## Outputs
main.log	- same output as that on the terminal while running the CPP program
polymer.XYZ 	- type and cartesian coordinates of all particles, can be opened in Ovito and XMD
polymer.XSF	- box size, and coordinates of all particles, particles are not distinguishable as Carbon is used to represent both the types, can be opened in VESTA
polymer.data	- data file in the prescribed LAMMPS format, can be read in using 'read_data' command in LAMMPS,

## TODO

- [ ] arbitrary constraints on a particular component of positions (eg. grow polymer backbone on a plane)
- [x] Use Classes/Objects to pass the data around
- [ ] Input file for reading the customizable parameters
- [ ] Periodic Lookup table
- [ ] Construct backbone and then add sideChains