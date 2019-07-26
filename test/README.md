1. Modify and execute using the input file
	`sarw -i pe.in | tee pe.out`
2. Open the output files using Ovito, VESTA or LAMMPS.
	For example,
	- `ovito polymer.dat` and select **molecular** as atomstyle
	- `ovito polymer.xyz`
	- `VESTA polymer.xsf`
	- in LAMMPS use `read_data polymer.dat`
