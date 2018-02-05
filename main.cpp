/*
*	Developed by Abhishek Shandilya
*	to-do
	1. periodic lookup
	2. book-keeping bonds
*/

// latest strategy - 29 Jan 2018 - construct the backbone and then add side-branches


#include "Random.h"
#include "matrix3.h"
#include "scalar.h"
#include "Euler.h"
#include "LatticeUtils.h"
#include <iostream>
#include <fstream>
#include "main.h"

const char species[6][10] = {"", "CH3", "CH2", "CH", "CH_aro", "C_aro"};

const double bondLengths[] =
	{
		1.54,	//CH3 	 - CH2 		or CH2 - CH
		1.51,	//CH 	 - C_aro
		1.40	//CH_aro - CH_aro
	};

// const double bondAngles[] =
// 	{
// 		112,	//CH3 -	CH - CH2
// 		114,	//CH - CH2 - CH
// 		109.5,	//CH3 - CH - C_aro
// 		120,	//C_aro - CH_aro - CH_aro
// 		120,	//CH_aro - CH_aro - CH_aro
// 		120,	//CH_aro - CH_aro - C_aro
// 		120		//CH_aro - C_aro - CH
// 	};

const double atomMass[] = {15, 14, 13, 13, 12};

struct unitedAtom
{
	int type;
	/*
		types:
		1. CH3
		2. CH2
		3. CH
	*/
	vector3<> pos;
	unitedAtom() : type(0), pos(vector3<>(0,0,0)) {}
	unitedAtom(int type, vector3<> pos) : type(type), pos(pos) {}
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

// std::vector<Bond> listBonds;
// std::vector<Angle> listAngles;
// std::vector<Dihedral> listDihedrals(nAtoms, Dihedral());
// std::vector<Dihedral> listImpropers(nAtoms, Dihedral());

// helper functions
vector3<> randomUnitStep();
void tetrahedral(const vector3<>, const vector3<>, vector3<>&, vector3<>&);
vector3<> wrapPBC(const vector3<>);
vector3<> randomConePos(const std::vector<unitedAtom>, const int, const int, const double);

bool checkCollision(const std::vector<unitedAtom>, const std::vector<unitedAtom>, const int ignoreIndex = -1);

// export functions
void exportXYZ(const std::vector<unitedAtom>);
void exportXSF(const std::vector<unitedAtom>);
void exportLAMMPS(const std::vector<unitedAtom>);
void printReport(const std::vector<unitedAtom>, std::vector<int>);
int minElement(const std::vector<int> v);
int maxElement(const std::vector<int> v);

// chain operations
bool initiateChain(std::vector<unitedAtom>&, std::vector<int>&, std::vector<int>&, std::vector<int> &);
// bool addSecondHalf(std::vector<unitedAtom>&, const std::vector<unitedAtom>);
bool propagateChain(std::vector<unitedAtom>&, std::vector<int>&, std::vector<int>&, std::vector<int> &);
void terminateChain(std::vector<unitedAtom>&);
bool randomSeed(std::vector<unitedAtom>&, const std::vector<unitedAtom>);

// void initialise(std::vector<unitedAtom> &);

int main(int argc, char *argv[])
{
	std::vector<unitedAtom> polymerChains;	// contains all finalised united atoms
	std::vector<int> lastIndices(nChains);
	std::vector<int> penultimateIndices(nChains);
	std::vector<int> chainLengths(nChains);
	
	printf("===Self Avoiding Random Walk===\n");
	
	int shortestChain = 0, longestChain = 0;
	if (initiateChain(polymerChains, lastIndices, penultimateIndices, chainLengths))
		while (polymerChains.size() < nUnitedAtoms)
		{
			if (!propagateChain(polymerChains, lastIndices, penultimateIndices, chainLengths)) break;
			// shortestChain = minElement(chainLengths);
			// longestChain = minElement(chainLengths);
			// if (DEBUG) printf("DEBUG:: longestChain: %d\n", longestChain);
		}
	terminateChain(polymerChains);

	printf("\n===Exporting Files===\n");
	// export to XYZ file for OVITO
	exportXYZ(polymerChains);
	// export to XSF file for VESTA
	exportXSF(polymerChains);
	// export to LAMMPS data file
	exportLAMMPS(polymerChains);

	printReport(polymerChains, chainLengths);

	return 0;
}

/* 
	initialise 
		- the side-branch based on desired polymer
	useless for now
*/
void initialise(std::vector<unitedAtom> &sideBranch)
{

	// Aromatic Ring for poly-styrene
	#ifdef POLYSTYRENE
	double d  = bondLengths[2];
	double dz = d * cos(M_PI/3);
	double dx = d * sin(M_PI/3);

	sideBranch.push_back( unitedAtom(5, vector3<>(  0, 0,    0)) );
	sideBranch.push_back( unitedAtom(4, vector3<>(-dx, 0,   dz)) );
	sideBranch.push_back( unitedAtom(4, vector3<>(-dx, 0,   dz + d)) );
	sideBranch.push_back( unitedAtom(4, vector3<>(  0, 0, 2*dz + d)) );
	sideBranch.push_back( unitedAtom(4, vector3<>( dx, 0,   dz + d)) );
	sideBranch.push_back( unitedAtom(4, vector3<>( dx, 0,   dz)) );
	#endif

	// CH3 for Propylene
	#ifdef POLYPROPYLENE
	sideBranch.push_back( unitedAtom(1, vector3<>(  0, 0,    0)) );
	#endif

	// nothing for Poly-ethylene
	#ifdef POLYETHYLENE
	// sideBranch.push_back( unitedAtom(1, vector3<>(  0, 0,    0)) );
	#endif
}

vector3<> randomUnitStep()
{
	vector3<> step(Random::uniform(), Random::uniform(), Random::uniform());
	step = normalize(step);
	return step;
}


int maxElement(const std::vector<int> v)
{
	int ans = v[0];
	for (int i : v)
		if (ans < i) ans = i;
	return ans;
}

int minElement(const std::vector<int> v)
{
	int ans = v[0];
	for (int i : v)
		if (ans > i) ans = i;
	return ans;
}

// useless for now
void tetrahedral(const vector3<> a, const vector3<> b, vector3<> &c, vector3<> &d)
{
	vector3<> z = normalize(cross(a, b));
	vector3<> y = normalize(a + b);
	vector3<> x = normalize(cross(y, z));

	matrix3<> rot1, rot2, rot_mat;

	rot1.set_rows(x, y, z);
	rot2 = matrixFromEuler(vector3<>(M_PI, M_PI/2, 0));
	rot_mat = ~rot1*rot2*rot1;

	c = rot_mat * a;
	d = rot_mat * b;
}

// useless for now
vector3<> wrapPBC(const vector3<> pos)
{
	vector3<> wrappedPos = pos;
	for (int i = 0; i < 3; ++i)
	{
		if 		(pos[i] < 0) 	   wrappedPos[i] += boxSize;
		else if (pos[i] > boxSize) wrappedPos[i] -= boxSize;
	}
	return wrappedPos;
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
		if (ignoreIndex > 0 && i==ignoreIndex) continue;
		for (int j = 0; j < newChainSize; ++j)
		{
			vector3<> dist = (polymerChains[i].pos - newChainLinks[j].pos);
			dist[0] -= boxSize * nearbyint(dist[0] / boxSize);
			dist[1] -= boxSize * nearbyint(dist[1] / boxSize);
			dist[2] -= boxSize * nearbyint(dist[2] / boxSize);
			if (dist.length() < minDist)
			{
				flag = true;
				break;
			}
		}
		if (flag) break;
	}
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
	fprintf(fp, "%d.0\t0.0\t0.0\n0.0\t%d.0\t0.0\n0.0\t0.0\t%d.0\n", boxSize, boxSize, boxSize);

	fprintf(fp, "CONVVEC\n");
	fprintf(fp, "%d.0\t0.0\t0.0\n0.0\t%d.0\t0.0\n0.0\t0.0\t%d.0\n", boxSize, boxSize, boxSize);

	fprintf(fp, "PRIMCOORD\n%d\t%d\n", polymerChains.size(), 1);
	for (unitedAtom atom : polymerChains)
		fprintf(fp, "C\t%lf\t%lf\t%lf\n", atom.pos[0], atom.pos[1], atom.pos[2]);
	fclose(fp);
}

void exportLAMMPS(const std::vector<unitedAtom> polymerChains)
{
	printf("Exporting LAMMPS data file: polymer.data\n");
	std::vector<Bond> 		listBonds;
	std::vector<Angle> 		listAngles;
	std::vector<Dihedral> 	listDihedrals;

	// std::vector<int> bondVector = {1, 2, 3, 3, 3, 3, 3, 3};
	// std::vector<int> angleVector = {109.5, 120, 120, 120};

	// ------------ Counting bonds ----------------
	for (int i = 0; i < polymerChains.size()-1; ++i)
	{
		if ((polymerChains[i].type == 1 && polymerChains[i+1].type == 1)) continue; // avoid the linking CH3 to another CH3
		listBonds.push_back(Bond(1, i+1, i+2));	// offset of 1 to start index at 1
	}
	// ------------ Counting angles ----------------
	for (int i = 0; i < polymerChains.size()-2; ++i)
	{
		if ((polymerChains[i  ].type == 1 && polymerChains[i+1].type == 1)) continue; 	// avoid the linking CH3 to another CH3
		if ((polymerChains[i+1].type == 1 && polymerChains[i+2].type == 1)) continue;
		listAngles.push_back(Angle(1, i+1, i+2, i+3));
	}
	// ------------ Counting dihedrals ----------------
	for (int i = 0; i < polymerChains.size()-3; ++i)
	{
		if ((polymerChains[i  ].type == 1 && polymerChains[i+1].type == 1)) continue; 	// avoid the linking CH3 to another CH3
		if ((polymerChains[i+1].type == 1 && polymerChains[i+2].type == 1)) continue;
		if ((polymerChains[i+2].type == 1 && polymerChains[i+3].type == 1)) continue;
		listDihedrals.push_back(Dihedral(1, i+1, i+2, i+3, i+4));
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
	fprintf(fp, "%lf %lf xlo xhi\n", -boxSize/2., boxSize/2.);
	fprintf(fp, "%lf %lf ylo yhi\n", -boxSize/2., boxSize/2.);
	fprintf(fp, "%lf %lf zlo zhi\n", -boxSize/2., boxSize/2.);
	fprintf(fp, "\n");

	fprintf(fp, "Masses\n\n");
	for (int i=0;i<2;i++) fprintf(fp, "%d\t%lf\n", i+1, atomMass[i]);
	fprintf(fp, "\n");
	// ========= Coefficients ==============

	// ========= Atoms ==============
	fprintf(fp, "Atoms\n\n");
	int i = 0;
	for (unitedAtom ua : polymerChains)
		fprintf(fp, "%d\t%d\t1\t%lf\t%lf\t%lf\n", ++i, ua.type, ua.pos[0], ua.pos[1], ua.pos[2]);
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
	
	printf("\nBoxSize\t: %d\n", boxSize);
	
	printf("\nnUnitedAtoms::\n");
	printf("Desired\t: %d\n", nUnitedAtoms);
	printf("Actual\t: %d\n", polymerChains.size());
	
	// printf("\nnChains::\n");
	// printf("Desired\t: %d\n", nChains);
	// printf("Actual\t: %d\n",  nChains_actual);
	
	printf("\nNumber Density::\n");
	printf("Desired\t: %.2f\n", pow(10,-24)*number_density);
	printf("Actual\t: %.2f\n", (polymerChains.size())/pow(boxSize, 3));

	printf("\nMass Density::\n");
	printf("Desired\t: %.2f\n", mass_density);
	printf("Actual\t: %.2f\n", ( 14 * polymerChains.size() )/(N_avogadro * pow(boxSize*pow(10,-8), 3)));

	// for (int i = 0; i < chainLengths.size(); i+=2)
	// 	chainLengths[i] = (chainLengths[i+1] - chainLengths[i]);

	printf("Exporting distribution of chain lengths: chains.dat\n");
	FILE* fp = fopen("chains.dat", "w");
	for (int i : chainLengths) fprintf(fp, "%d\n", i);
	fclose(fp);
}

// ring and CH2 together - useless for now
bool addSecondHalf(std::vector<unitedAtom> &polymerChains, const std::vector<unitedAtom> sideBranch)
{
	// initialize local variables
	vector3<> ringPos, r12, r10, rCH, rCH2, pos1, pos2;
	matrix3<> localCoords;
	int iTrial, index;

	index = polymerChains.size();
	rCH  = polymerChains[index-1].pos;
	
	// try adding aromatic ring and CH2
	std::vector<unitedAtom> newChainLinks;
	iTrial = 0;
	while (iTrial++ < maxTrials)
	{
		newChainLinks.clear();
		ringPos = randomConePos(polymerChains, index-1, index-2, bondLengths[1]);

		// calculate the plane to place the ring 
		r12 = normalize(polymerChains[index-2].pos - rCH);
		r10 = normalize(ringPos - rCH);
		localCoords.set_row(2, r10);
		localCoords.set_row(0, normalize(r12 - r10 * dot(r10, r12)));
		localCoords.set_row(1, normalize(cross(localCoords.row(2), localCoords.row(0))));

		for (int i = 0; i < sideBranch.size(); ++i)
			newChainLinks.push_back( unitedAtom(sideBranch[i].type, ringPos + sideBranch[i].pos * localCoords) );

		tetrahedral(r12, r10, pos1, pos2);
		// discrete choice
		rCH2 = (Random::uniform() > 0.5)? pos1 : pos2;
		newChainLinks.push_back( unitedAtom(2, bondLengths[1]*rCH2 + rCH) );

		if (checkCollision(polymerChains, newChainLinks)) continue;
		polymerChains.insert(polymerChains.end(), newChainLinks.begin(), newChainLinks.end());
		break;
	}

	if (iTrial < maxTrials) return true;
	
	iTrial = 0;
	while (iTrial++ < maxTrials)
	{
		newChainLinks.clear();
		ringPos = randomConePos(polymerChains, index-1, index-2, bondLengths[1]);

		// calculate the plane to place the ring 
		r12 = normalize(polymerChains[index-2].pos - rCH);
		r10 = normalize(ringPos - rCH);
		localCoords.set_row(2, r10);
		localCoords.set_row(0, normalize(r12 - r10 * dot(r10, r12)));
		localCoords.set_row(1, normalize(cross(localCoords.row(2), localCoords.row(0))));

		for (int i = 0; i < sideBranch.size(); ++i)
			newChainLinks.push_back( unitedAtom(sideBranch[i].type, ringPos + sideBranch[i].pos * localCoords) );

		if (checkCollision(polymerChains, newChainLinks)) continue;
		polymerChains.insert(polymerChains.end(), newChainLinks.begin(), newChainLinks.end());
		break;
	}
	// printf("Tried %d times for just the ring at %d\n", iTrial, index);
	// if (iTrial >= maxTrials) return false;
	
	// iTrial = 0;
	// newChainLinks.resize(1);
	// while (iTrial++ < maxTrials)
	// {
	// 	ringPos = randomConePos(polymerChains, index-1, index-2, bondLengths[1]);
	// 	r12 = normalize(polymerChains[index-2].pos - rCH);
	// 	r10 = normalize(ringPos - rCH);

	// 	tetrahedral(r12, r10, pos1, pos2);
	// 	// discrete choice
	// 	rCH2 = (Random::uniform() > 0.5)? pos1 : pos2;
	// 	newChainLinks[0] = unitedAtom(2, bondLengths[1]*rCH2 + rCH);

	// 	if (checkCollision(polymerChains, newChainLinks)) continue;
	// 	polymerChains.insert(polymerChains.end(), newChainLinks.begin(), newChainLinks.end());
	// 	break;
	// }
	// printf("Tried %d times for CH2 after the ring at %d\n", iTrial, index);
	// if (iTrial < maxTrials) return true;
	return false;
}

// add CH2 on a cone of tetrahedral angle
bool propagateChain(std::vector<unitedAtom> &polymerChains, std::vector<int> &lastIndices, std::vector<int> &penultimateIndices, std::vector<int> &chainLengths)
{
	if (DEBUG) printf("DEBUG:: Entered Propagation step\n");
	// propagate all chains, if possible
	// udpate indices and lengths
	int lastIndex, penultimateIndex;

	std::vector<unitedAtom> newCH2(1, unitedAtom());
	
	for (int i = 0; i < nChains; ++i)
	{
		lastIndex 		 = 		  lastIndices[i];
		penultimateIndex = penultimateIndices[i];
		int iTrial = 0;
		while (iTrial++ < maxTrials)
		{
			// avoid hard-coding any number
			newCH2[0] = unitedAtom(2, randomConePos(polymerChains, lastIndex, penultimateIndex, bondLengths[0]));
			if (checkCollision(polymerChains, newCH2, lastIndex)) continue;
			polymerChains.push_back(newCH2[0]);

			penultimateIndices[i] = lastIndices[i];
			lastIndices[i] 		  = polymerChains.size()-1;

			chainLengths[i]++;
			break;
		}
		if (DEBUG) printf("DEBUG:: iTrial %d\n", iTrial);
		if (DEBUG && iTrial-1 < maxTrials)
		{
			printf("DEBUG:: Propagated %dth chain at %d\n", i, polymerChains.size());
			// printf("%.3f\n", newCH2[0].pos[0]);
		}
	}
	// if (iTrial >= maxTrials){
	//	if (DEBUG) printf("DEBUG:: Propagatation failed at %d\n", polymerChains.size());
	//	return false;
	//}
	return true;
}

/*
Try adding a pair of CH3 and CH2 to initiate a new chain
	if successful - add the united-atoms to data-structure, and return true
	if failed - return false, which means no more chains can be added
*/
bool initiateChain(std::vector<unitedAtom> &polymerChains, std::vector<int> &lastIndices, std::vector<int> &penultimateIndices, std::vector<int> &chainLengths)
{
	bool flag = true;
	if (DEBUG) printf("DEBUG:: Entered Initiation step\n");
	std::vector<unitedAtom> newChainLinks(2, unitedAtom());
	for (int i = 0; i < nChains; ++i)
	{
		if (randomSeed(newChainLinks, polymerChains))
		{
			polymerChains.push_back(newChainLinks[0]);
			polymerChains.push_back(newChainLinks[1]);
			// update indices
			lastIndices[i] = polymerChains.size()-1;
			penultimateIndices[i] = polymerChains.size()-2;
			chainLengths[i]+=2;
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
			Random::uniform(0, boxSize),
			Random::uniform(0, boxSize),
			Random::uniform(0, boxSize));

		// CH atom at a random position on a sphere around the CH3 seed
		step = bondLengths[0] * randomUnitStep();
		newChainLinks[0] = (unitedAtom(1, pos0));
		newChainLinks[1] = (unitedAtom(2, pos0 + step));
		if (!checkCollision(polymerChains, newChainLinks)) break;
	}
	if (iTrial >= maxTrials) return false;
	return true;
}

void terminateChain(std::vector<unitedAtom> &polymerChains)
{
	int chainSize = polymerChains.size();

	if (DEBUG) printf("DEBUG:: Terminated!! at %d\n", chainSize);
	switch(polymerChains[chainSize-1].type)
	{
		case 2:
			polymerChains[chainSize-1].type = 1;
			break;
		// case 5: // can't add the CH2
		// 	break;
		// case 3:	// can't add the ring and/or CH2
		// 	if (chainSize >= 2)
		// 	{
		// 		polymerChains.pop_back();
		// 		polymerChains.pop_back();
		// 		polymerChains[chainSize-9].type = 2;
		// 	}
		// 	break;
	}
}
