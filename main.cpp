/*
*	Developed by Abhishek Shandilya
*	to-do
	1. periodic lookup
	2. book-keeping bonds
*/

// latest strategy - 29 Jan 2018 - construct the backbone and then add side-branches

#include <core/Random.h>
#include <core/matrix3.h>
#include <fluid/Euler.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "main.h"

const char species[6][10] = {"", "CH3", "CH2", "CH", "CH_aro", "C_aro"};

const double bondLengths[] =
	{
		1.54,	//CH3 	 - CH2 		or CH2 - CH
		1.51,	//CH 	 - C_aro
		1.40	//CH_aro - CH_aro
	};

const double atomMass[] = {15, 14, 13, 13, 12};

struct unitedAtom
{
	int type;
	/*
		types:
		1. CH3
		2. CH2
		3. CH
		4. Si
	*/
	int chainID;
	vector3<> pos;
	unitedAtom() : type(0), chainID(0), pos(vector3<>(0,0,0)) {}
	unitedAtom(int type, vector3<> pos) : type(type), pos(pos) {}
	unitedAtom(int type, int chainID, vector3<> pos) : type(type), chainID(chainID), pos(pos) {}
};

struct Bond
{
	int type;
	/*
		types:
		1. CH3-CH2 or CH2-CH2
	*/
	int species[2];
	Bond(int t=0, int a=0, int b=0)
	{
		type = t;
		species[0] = a; species[1] = b;
	}
};

struct Angle
{
	int type;
	// 1. 109.5
	// 2. 120
	int species[3];
	Angle(int t=0, int a=0, int b=0, int c=0)
	{
		type = t;
		species[0] = a; species[1] = b; species[2] = c;
	}
};

struct Dihedral
{
	int type;
	int species[4];
	Dihedral(int t=0, int a=0, int b=0, int c=0, int d=0)
	{
		type = t;
		species[0] = a; species[1] = b; species[2] = c; species[3] = d;
	}
};

// helper functions
vector3<> randomUnitStep();
vector3<> randomConePos(const std::vector<unitedAtom>, const int, const int, const double);

bool checkCollision(const std::vector<unitedAtom>, const std::vector<unitedAtom>, const int ignoreIndex = -1);

// export functions
void exportXYZ(const std::vector<unitedAtom>);
void exportXSF(const std::vector<unitedAtom>);
void exportLAMMPS(const std::vector<unitedAtom>);
void printReport(const std::vector<unitedAtom>, std::vector<int>);

// chain operations
bool initiateChain(std::vector<unitedAtom>&, std::vector<int>&, std::vector<int>&, std::vector<int> &);
bool propagateChain(std::vector<unitedAtom>&, std::vector<int>&, std::vector<int>&, std::vector<int> &);
void terminateChain(std::vector<unitedAtom>&, const std::vector<int>);
bool randomSeed(std::vector<unitedAtom>&, const std::vector<unitedAtom>);

int main(int argc, char *argv[])
{
	std::vector<unitedAtom> polymerChains;	// contains all finalised united atoms
	std::vector<int> lastIndices(nChains); // index of last united atom of each chain in polymerChains[]
	std::vector<int> penultimateIndices(nChains); // index of penulutimate united atom of each chain in polymerChains[]
	std::vector<int> chainLengths(nChains);	// length of each polymer chain

	#ifdef CENTER_SEED
	int seed_vol = (4./3) * M_PI * pow(seed_radius, 3);
	boxSize = nearbyint(pow(pow(boxSize, 3) + seed_vol, 0.3333));
	unitedAtom center_seed = unitedAtom(4, vector3<>(boxSize*0.5, boxSize*0.5, boxSize*0.5));
	// polymerChains.push_back(center_seed);
	#endif
	
	printf("===Self Avoiding Random Walk===\n");
	
	int shortestChain = 0, longestChain = 0;
	if (initiateChain(polymerChains, lastIndices, penultimateIndices, chainLengths))
		while (polymerChains.size() < nUnitedAtoms)
			if (!propagateChain(polymerChains, lastIndices, penultimateIndices, chainLengths)) break;
	terminateChain(polymerChains, lastIndices);

	printf("\n===Exporting Files===\n");
	exportXYZ(polymerChains); // export to XYZ file for OVITO
	exportXSF(polymerChains); // export to XSF file for VESTA
	exportLAMMPS(polymerChains); // export to LAMMPS data file

	printReport(polymerChains, chainLengths);

	return 0;
}

vector3<> randomUnitStep()
{
	vector3<> step(Random::uniform(), Random::uniform(), Random::uniform());
	step = normalize(step);
	return step;
}

vector3<> randomConePos(const std::vector<unitedAtom> polymerChains, const int lastIndex, const int penultimateIndex, const double distance)
{
	double theta, phi;
	vector3<> localX, localY, localZ, newPos;

	phi = (180 - 109.5)*M_PI/180;

	// calculate local coordinate system for the cone
	localZ = normalize(polymerChains[lastIndex].pos - polymerChains[penultimateIndex].pos);

	if (abs(localZ[0]) > 0) localX = vector3<>(1,0,0) - localZ * localZ.x();
	else					localX = vector3<>(0,1,0) - localZ * localZ.y();
	localX = normalize(localX);
	localY = normalize(cross(localZ, localX));

	theta = Random::uniform(0, 359)*M_PI/180;

	newPos = distance*sin(phi)*cos(theta)*localX
	       + distance*sin(phi)*sin(theta)*localY
		   + distance*cos(phi)		     *localZ;
	newPos += polymerChains[lastIndex].pos;
	return newPos;
}

bool checkCollision(
	const std::vector<unitedAtom> polymerChains,
	const std::vector<unitedAtom> newChainLinks,
	const int ignoreIndex)
{
	bool flag = false;
	int iMax, newChainSize, chainSize;

	chainSize 	 = polymerChains.size();
	newChainSize = newChainLinks.size();
	
	for (int i = 0; i < chainSize; ++i)
	{
		if (ignoreIndex >= 0 && i==ignoreIndex) continue;
		for (int j = 0; j < newChainSize; ++j)
		{
			vector3<> dist = (polymerChains[i].pos - newChainLinks[j].pos);
			for (int k = 0; k < 3; ++k)
				dist[k] -= boxSize[k] * nearbyint(dist[k] / boxSize[k]);
			if (dist.length() < minDist)
			{
				flag = true;
				break;
			}
		}
		if (flag) break;
	}
	if (flag) return flag;

	#ifdef CENTER_SEED
	for (int j = 0; j < newChainSize; ++j)
	{
		unitedAtom center_seed = unitedAtom(4, vector3<>(boxSize[0]*0.5, boxSize[1]*0.5, boxSize[2]*0.5));
		vector3<> dist = (newChainLinks[j].pos - center_seed.pos);
		if (dist.length() < seed_radius) { flag = true; break;}
	}
	#endif

	return flag;
}

void exportXYZ(const std::vector<unitedAtom> polymerChains)
{
	printf("Exporting XYZ file: polymer.xyz\n");
	FILE* fp = fopen("polymer.xyz", "w");
	fprintf(fp, "%d\n", polymerChains.size());
	fprintf(fp, "POLYSTYRENE\n");
	for (unitedAtom atom : polymerChains)
		fprintf(fp, "%d\t%lf\t%lf\t%lf\n", atom.type, atom.pos[0], atom.pos[1], atom.pos[2]);
	fclose(fp);
}

void exportXSF(const std::vector<unitedAtom> polymerChains)
{
	printf("Exporting XSF file: polymer.xsf\n");
	FILE* fp = fopen("polymer.xsf", "w");
	fprintf(fp, "CRYSTAL\nPRIMVEC\n");
	fprintf(fp, "%d.0\t0.0\t0.0\n0.0\t%d.0\t0.0\n0.0\t0.0\t%d.0\n", boxSize[0], boxSize[1], boxSize[2]);

	fprintf(fp, "CONVVEC\n");
	fprintf(fp, "%d.0\t0.0\t0.0\n0.0\t%d.0\t0.0\n0.0\t0.0\t%d.0\n", boxSize[0], boxSize[1], boxSize[2]);

	fprintf(fp, "PRIMCOORD\n%d\t%d\n", polymerChains.size(), 1);
	for (unitedAtom atom : polymerChains)
		fprintf(fp, "C\t%lf\t%lf\t%lf\n", atom.pos[0], atom.pos[1], atom.pos[2]);
	fclose(fp);
}

std::vector<Bond> 		listBonds;
void exportLAMMPS(const std::vector<unitedAtom> polymerChains)
{
	printf("Exporting LAMMPS data file: polymer.data\n");
	std::vector<Angle> 		listAngles;
	std::vector<Dihedral> 	listDihedrals;

	// ------------ Counting bonds ----------------
	// ------------ Counting angles ----------------
	int a1, a2, a3, a4;
	int i = 0;
	for (int iChain = 1; iChain <= nChains; ++iChain)
	{
		i = 2*iChain - 2;
		// find the very first monomer
		a1 = i;
		// find the second monomer
		while(++i < polymerChains.size())
		{
			if (polymerChains[i].chainID == iChain)
			{
				a2 = i;
				break;
			}
		}
		// find the third monomer
		while(++i < polymerChains.size())
		{
			if (polymerChains[i].chainID == iChain)
			{
				a3 = i;
				break;
			}
		}
		while(i < polymerChains.size())
		{
			// add to angles list
			listAngles.push_back(Angle(1, a1+1, a2+1, a3+1));
			// udpate all three and repeat
			a1 = a2; a2 = a3;
			while(++i < polymerChains.size())
			{
				if (polymerChains[i].chainID == iChain)
				{
					a3 = i;
					break;
				}
			}
		}
	}
	// ------------ Counting dihedrals ----------------
	for (int iChain = 1; iChain <= nChains; ++iChain)
	{
		i = 2*iChain - 2;
		// find the very first monomer
		a1 = i;
		// find the second monomer
		while(++i < polymerChains.size())
		{
			if (polymerChains[i].chainID == iChain)
			{
				a2 = i;
				break;
			}
		}
		// find the third monomer
		while(++i < polymerChains.size())
		{
			if (polymerChains[i].chainID == iChain)
			{
				a3 = i;
				break;
			}
		}
		// find the fourth monomer
		while(++i < polymerChains.size())
		{
			if (polymerChains[i].chainID == iChain)
			{
				a4 = i;
				break;
			}
		}
		while(i < polymerChains.size())
		{
			// add to dihedrals list
			listDihedrals.push_back(Dihedral(1, a1+1, a2+1, a3+1, a4+1));
			// udpate all four and repeat
			a1 = a2; a2 = a3; a3 = a4;
			while(++i < polymerChains.size())
			{
				if (polymerChains[i].chainID == iChain)
				{
					a4 = i;
					break;
				}
			}
		}
	}
	// ------------ Counting impropers ----------------

	FILE* fp = fopen("polymer.data", "w");
	#ifdef POLYSTYRENE
	fprintf(fp, "POLYSTYRENE CHAINS\n\n");
	#endif
	#ifdef POLYPROPYLENE
	fprintf(fp, "POLYPROPYLENE CHAINS\n\n");
	#endif
	#ifdef POLYETHYLENE
	fprintf(fp, "POLYETHYLENE CHAINS\n\n");
	#endif

	fprintf(fp, "%5d\tatoms\n", polymerChains.size());
	fprintf(fp, "%5d\tbonds\n", listBonds.size());
	fprintf(fp, "%5d\tangles\n", listAngles.size());
	fprintf(fp, "%5d\tdihedrals\n", listDihedrals.size());
	// fprintf(fp, "%d\timpropers\n", nImpropers);
	fprintf(fp, "\n");
	fprintf(fp, "2\tatom types\n");
	fprintf(fp, "1\tbond types\n");
	fprintf(fp, "1\tangle types\n");
	fprintf(fp, "1\tdihedral types\n");
	// fprintf(fp, "1\timproper types\n");
	fprintf(fp, "\n");
	fprintf(fp, "0.0 %d xlo xhi\n", boxSize[0]);
	fprintf(fp, "0.0 %d ylo yhi\n", boxSize[1]);
	fprintf(fp, "0.0 %d zlo zhi\n", boxSize[2]);
	fprintf(fp, "\n");

	fprintf(fp, "Masses\n\n");
	for (int i=0;i<2;i++) fprintf(fp, "%d\t%lf\n", i+1, atomMass[i]);
	fprintf(fp, "\n");

	// ========= Atoms ==============
	fprintf(fp, "Atoms\n\n");
	i = 0;
	for (unitedAtom ua : polymerChains)
		fprintf(fp, "%d\t%d\t%d\t%lf\t%lf\t%lf\n", ++i, ua.chainID, ua.type, ua.pos[0], ua.pos[1], ua.pos[2]);
	fprintf(fp, "\n");

	fprintf(fp, "Bonds\n\n");
	i = 0;
	for(Bond b : listBonds)
		fprintf(fp, "%d\t%d\t%d\t%d\n", ++i, b.type, b.species[0], b.species[1]);
	fprintf(fp, "\n");
	
	fprintf(fp, "Angles\n\n");
	i = 0;
	for(Angle b : listAngles)
		fprintf(fp, "%d\t%d\t%d\t%d\t%d\n", ++i, b.type, b.species[0], b.species[1], b.species[2]);
	fprintf(fp, "\n");
	
	fprintf(fp, "Dihedrals\n\n");
	for (int i = 0; i < listDihedrals.size(); i++)
	{
		fprintf(fp, "%d\t%d\t", i+1, listDihedrals[i].type);
		for (int j = 0; j < 4; ++j)	fprintf(fp, "%d\t", listDihedrals[i].species[j]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	
	// fprintf(fp, "Impropers\n");
	// for (int i = 0; i < nImpropers; i++)
	// {
	// 	fprintf(fp, "%d\t%d\t", i+1, listImpropers[i].type);
	// 	for (int j = 0; j < 4; ++j)	fprintf(fp, "%d\t", listImpropers[i].species[j]);
	// 	fprintf(fp, "\n");
	// }
	// fprintf(fp, "\n");
	fclose(fp);
}

void printReport(const std::vector<unitedAtom> polymerChains, const std::vector<int> chainLengths)
{
	int iAtom = 0, nChains_actual = 0;
	// std::vector<int> chainLengths;
	
	// for(unitedAtom i : polymerChains)
	// {
	// 	iAtom++;
	// 	if (i.type == 1)
	// 	{
	// 		nChains_actual++;
	// 		chainLengths.push_back(iAtom);
	// 	}
	// }
	// nChains_actual /= 2;

	printf("\n===Report===\n");
	
	printf("\nBoxSize\t: %d x %d x %d\n", boxSize[0], boxSize[1], boxSize[2]);
	
	printf("\nnUnitedAtoms::\n");
	printf("Desired\t: %d\n", nUnitedAtoms);
	printf("Actual\t: %d\n", polymerChains.size());
	
	// printf("\nnChains::\n");
	// printf("Desired\t: %d\n", nChains);
	// printf("Actual\t: %d\n",  nChains_actual);
	
	int vol = (boxSize[0]*boxSize[1]*boxSize[2]);
	printf("\nNumber Density::\n");
	printf("Desired\t: %.3e\n", number_density);
	printf("Actual\t: %.3e\n", polymerChains.size()/(vol*pow(10,-24)));

	printf("\nMass Density::\n");
	printf("Desired\t: %.2f\n", mass_density);
	printf("Actual\t: %.2f\n", ( 14*polymerChains.size() )/(N_avogadro*vol*pow(10,-24)) );

	// for (int i = 0; i < chainLengths.size(); i+=2)
	// 	chainLengths[i] = (chainLengths[i+1] - chainLengths[i]);

	printf("Exporting distribution of chain lengths: chains.dat\n");
	FILE* fp = fopen("chains.dat", "w");
	for (int i : chainLengths) fprintf(fp, "%d\n", i);
	fclose(fp);
}

// add CH2 on a cone of tetrahedral angle
bool propagateChain(std::vector<unitedAtom> &polymerChains, std::vector<int> &lastIndices, std::vector<int> &penultimateIndices, std::vector<int> &chainLengths)
{
	if (DEBUG) printf("DEBUG:: Entered Propagation step\n");
	// propagate all chains, if possible
	// udpate indices and lengths
	int iChain, lastIndex, penultimateIndex;
	bool flag = true, lengthFlag;

	std::vector<unitedAtom> newCH2(1, unitedAtom());
	
	for (int i = 0; i < nChains; ++i)
	{
		iChain = nearbyint(Random::uniform(0, nChains-1));
		lastIndex 		 = 		  lastIndices[iChain];
		penultimateIndex = penultimateIndices[iChain];
		int iTrial = 0;
		while (iTrial++ < maxTrials)
		{
			// avoid hard-coding any number
			newCH2[0] = unitedAtom(2, iChain+1, randomConePos(polymerChains, lastIndex, penultimateIndex, bondLengths[0]));
			if (checkCollision(polymerChains, newCH2, lastIndex)) continue;
			lengthFlag = false;
			// lengthFlag = chainLengths[iChain] >= minChainLength;
			if (!lengthFlag)
			{
				polymerChains.push_back(newCH2[0]);

				penultimateIndices[iChain] = lastIndices[iChain];
				lastIndices[iChain] 	   = polymerChains.size() - 1;

				listBonds.push_back(Bond(1, penultimateIndices[iChain] + 1, lastIndices[iChain] + 1));

				chainLengths[iChain]++;
			}
			break;
		}
		flag = flag && (iTrial-1 >= maxTrials);
		// flag = flag && ((iTrial-1 >= maxTrials) || lengthFlag);

		if (DEBUG) printf("DEBUG:: iTrial %d\n", iTrial);
		if (iTrial-1 < maxTrials)
		{
			if (DEBUG) printf("DEBUG:: Propagated %dth chain at %d\n", i, polymerChains.size());
			if (LOG) printf("LOG:: Progress = %d\n", 100*polymerChains.size()/nUnitedAtoms);
		}
	}
	// if (iTrial >= maxTrials){
	//	if (DEBUG) printf("DEBUG:: Propagatation failed at %d\n", polymerChains.size());
	//	return false;
	//}
	if (flag) return false;
	return true;
}

/*
Try adding a pair of CH3 and CH2 to initiate a new chain
	if successful - add the united-atoms to data-structure, and return true
	if failed - return false, which means no more chains can be added
*/
bool initiateChain(std::vector<unitedAtom> &polymerChains, std::vector<int> &lastIndices, std::vector<int> &penultimateIndices, std::vector<int> &chainLengths)
{
	if (DEBUG) printf("DEBUG:: Entered Initiation step\n");

	bool flag = true;
	int index;
	std::vector<unitedAtom> newChainLinks(2, unitedAtom());

	for (int i = 0; i < nChains; ++i)
	{
		if (randomSeed(newChainLinks, polymerChains))
		{
			newChainLinks[0].chainID = i+1;
			newChainLinks[1].chainID = i+1;
			polymerChains.push_back(newChainLinks[0]);
			polymerChains.push_back(newChainLinks[1]);
			// update indices, index starts from 0
			lastIndices[i] = polymerChains.size() - 1;
			penultimateIndices[i] = polymerChains.size() - 2;
			// update listBond, index starts from 1
			listBonds.push_back(Bond(1, penultimateIndices[i]+1, lastIndices[i]+1));
			chainLengths[i] += 2;
			if (DEBUG) printf("DEBUG:: Initiated %dth seed\n", i);
		}
		else
		{
			flag = false;
			break;
		}
	}
	if (DEBUG && flag) printf("DEBUG:: Initiated %d seeds\n", nChains);

	return flag;
}

bool randomSeed(std::vector<unitedAtom> &newChainLinks, const std::vector<unitedAtom> polymerChains)
{
	vector3<> step, pos0;
	int iTrial = 0;
	while (iTrial++ < maxTrials)
	{
		// random CH3 seed in the box
		pos0 = vector3<>(
			Random::uniform(0, boxSize[0]),
			Random::uniform(0, boxSize[1]),
			0.0);
		if (RANDOM_SEED) pos0[2] = Random::uniform(0, boxSize[2]);

		// CH atom at a random position on a sphere around the CH3 seed
		step = bondLengths[0] * randomUnitStep();
		newChainLinks[0] = (unitedAtom(1, pos0));
		newChainLinks[1] = (unitedAtom(2, pos0 + step));
		if (!checkCollision(polymerChains, newChainLinks)) break;
	}
	if (iTrial >= maxTrials) return false;
	return true;
}

void terminateChain(std::vector<unitedAtom> &polymerChains, const std::vector<int> lastIndices)
{
	int chainSize = polymerChains.size();

	if (DEBUG) printf("DEBUG:: Terminated!! at %d\n", chainSize);
	for (int i = 0; i < nChains; ++i)
	{
		if (polymerChains[lastIndices[i]].type == 2)
		{
			polymerChains[lastIndices[i]].type = 1;
			if (DEBUG) printf("Terminated chain # %d\n", i+1);
		}
	}
}
