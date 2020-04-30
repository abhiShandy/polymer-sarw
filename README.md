# Overview

Using Self-Avoiding Random Walk (SARW) to generate polymer chains for MD simulations.
Each monomer is modelled using united-atoms, where all hydrogens are contracted into the carbon atoms.  For example, CH3 and CH2 are modelled as a single particle of different kinds, CH3 is a particle of mass 15 a.u. and CH2 has mass of 14 a.u.
Each subsequent united atom is placed randomly on a cone of half-angle = (180 - tetrahedral angle)=75.5 degress to capture the sp3 bond angle.

## Description of Model Parameters

All model parameters are in the sarw.h file

<table>
	<thead>
		<th>Parameter</th>
		<th>Description</th>
		<th>Possible values</th>
		<th>Default</th>
	</thead>
	<tr>
		<td>boundary</td>
		<td>Turn off periodic boundary conditions on one of the cartesian directions (similar to LAMMPS syntax).</td>
		<td>3 options: ppf, pfp, fpp</td>
		<td>ppp</td>
	</tr>
	<tr>
		<td>boxSize</td>
		<td>Size of the cuboid shaped simulation box in Angstroms</td>
		<td>Any set of 3 integers separated by commas</td>
		<td>10,10,10</td>
	</tr>
	<tr>
		<td>graftedSeeds</td>
		<td>If the grafted seeds need to be on a specific plane, or specific coordinates, then provide the appropriate argument. Example: 'x' means seeds are placed at `x=0` plane. 'file' means seed locations should be read from 'graft.dat' in the same directory as input file. By default, grafted seeds are randomly placed in the simulation box which should be used when no chains need to be grafted.</td>
		<td>x, y, z, zz, file, random</td>
		<td>random</td>
	</tr>
	<tr>
		<td>graftFraction</td>
		<td>Fraction of total chains which need to be grafted at locations specified by the `graftedSeeds` option.</td>
		<td>Any fraction less than 1</td>
		<td>0.0</td>
	</tr>
	<tr>
		<td>growthBias</td>
		<td>A vector used to grow the chains in specific cartesian direction.</td>
		<td>Any set of 3 integers separated by commas. Each bias can either be 0 +1 or -1</td>
		<td>0,0,0</td>
	</tr>
	<tr>
		<td>logProgress</td>
		<td>Whether to show overall progress of chain growth in the output file.</td>
		<td>yes/no</td>
		<td>yes</td>
	</tr>
	<tr>
		<td>logSteps</td>
		<td>Whether to show details about every Monte Carlo attempt, and each function call to grow the chains.</td>
		<td>yes/no</td>
		<td>no</td>
	</tr>
	<tr>
		<td>maxAtoms</td>
		<td>Set the upper cap on number of atoms in the simulation box. If set to -1, it is fixed based on `targetDensity`</td>
		<td>Any integer</td>
		<td>-1</td>
	</tr>
	<tr>
		<td>maxTrials</td>
		<td>Set an upper limit on the Monte Carlo trials.</td>
		<td>Any integer</td>
		<td>100</td>
	</tr>
	<tr>
		<td>minDist</td>
		<td>Set a minimum distance allowed between non-bonded united-atoms.</td>
		<td>Any positive floating number</td>
		<td>2.0</td>
	</tr>
	<tr>
		<td>nChains</td>
		<td>Set the number of chains that are grown in the simulation box. Some of them maybe grafted if `graftFraction` is non-zero.</td>
		<td>Any integer</td>
		<td>1</td>
	</tr>
	<tr>
		<td>polymer</td>
		<td>Name of the polymer used in the header of LAMMPS data file</td>
		<td>Any string</td>
		<td>polymer</td>
	</tr>
	<tr>
		<td>roundRobin</td>
		<td>Whether the chains should be grown in round-robin fashion or randomly</td>
		<td>yes/no</td>
		<td>yes</td>
	</tr>
	<tr>
		<td>targetMassDensity</td>
		<td>Set the mass density of simulation box in gram per cubic centimeter</td>
		<td>Any positive fraction</td>
		<td>0.92</td>
	</tr>

</table>

## Outputs
main.log	- same output as that on the terminal while running the CPP program
polymer.XYZ 	- type and cartesian coordinates of all particles, can be opened in Ovito and XMD
polymer.XSF	- box size, and coordinates of all particles, particles are not distinguishable as Carbon is used to represent both the types, can be opened in VESTA
polymer.data	- data file in the prescribed LAMMPS format, can be read in using 'read_data' command in LAMMPS,

## TODO

- [ ] arbitrary constraints on a particular component of positions (eg. grow polymer backbone on a plane)
- [x] Use Classes/Objects to pass the data around
- [x] Input file for reading the customizable parameters
- [x] Periodic Lookup table
- [ ] Construct backbone and then add sideChains
- [ ] Total time taken
