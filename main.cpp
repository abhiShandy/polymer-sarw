/*
*	Developed by Abhishek Shandilya
*/


#include "Random.h"
#include "matrix3.h"
#include "scalar.h"
#include "Euler.h"
#include "LatticeUtils.h"
#include <iostream>
#include <fstream>

const char species[6][10] = {"", "CH3", "CH2", "CH", "CH_aro", "C_aro"};

const double bondLengths[] =
	{
		1.54,	//CH3 	 - CH2 		or CH2 - CH
		1.51,	//CH 	 - C_aro
		1.40	//CH_aro - CH_aro
	};

const double bondAngles[] =
	{
		112,	//CH3 -	CH - CH2
		114,	//CH - CH2 - CH
		109.5,	//CH3 - CH - C_aro
		120,	//C_aro - CH_aro - CH_aro
		120,	//CH_aro - CH_aro - CH_aro
		120,	//CH_aro - CH_aro - C_aro
		120		//CH_aro - C_aro - CH
	};

const double atomMass[] = {15, 14, 13, 13, 12};

struct unitedAtom
{
	int type;
	vector3<> pos;
	unitedAtom() : type(0), pos(vector3<>(0,0,0)) {}
	unitedAtom(int type, vector3<> pos) : type(type), pos(pos) {}
};

struct Bond
{
	int type;
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

const int nMonomers = 50;	// number of monomers in the simulation box, calculate using desired density
#define POLYSTYRENE
// #define POLYPROPYLENE
std::vector<int> chainLengths;

// std::vector<Bond> listBonds;
// std::vector<Angle> listAngles;
// std::vector<Dihedral> listDihedrals(nAtoms, Dihedral());
// std::vector<Dihedral> listImpropers(nAtoms, Dihedral());

const double minDist = 2;
const int maxTrials = 100;
const vector3<> boxSize(10, 10, 10);

// helper functions
vector3<> randomUnitStep();
void tetrahedral(const vector3<>, const vector3<>, vector3<>&, vector3<>&);
vector3<> wrapPBC(vector3<>);
vector3<> randomConePos(const std::vector<unitedAtom>, const int, const int, const double);

bool checkCollision(const std::vector<unitedAtom>, const std::vector<unitedAtom>, const bool flag = false);

// export functions
void exportXYZ(const std::vector<unitedAtom>);
void exportXSF(const std::vector<unitedAtom>);
// void exportLAMMPS(const std::vector<unitedAtom>&);
void exportChainLengths();

// chain operations
bool initiateChain(std::vector<unitedAtom>&);
bool addSecondHalf(std::vector<unitedAtom>&, const std::vector<unitedAtom>);
bool addFirstHalf(std::vector<unitedAtom>&);
void terminateChain(std::vector<unitedAtom>&);

void initialise(std::vector<unitedAtom> &);

int main(int argc, char *argv[])
{
	std::vector<unitedAtom> polymerChains;
	std::vector<unitedAtom> sideBranch;

	printf("===Self Avoiding Random Walk===\n");
	initialise(sideBranch);
	
	int iMonomer = 0;

	if(initiateChain(polymerChains))
		while(iMonomer++ < nMonomers)
		{
			if(addSecondHalf(polymerChains, sideBranch))
				if(addFirstHalf(polymerChains)) continue;
			terminateChain(polymerChains);

			if(initiateChain(polymerChains)) continue;
			break;
		}
	terminateChain(polymerChains);

	// export to XYZ file
	exportXYZ(polymerChains);
	exportXSF(polymerChains);
	// export to LAMMPS data file
	// exportLAMMPS(polymerChains);

	exportChainLengths();

	return 0;
}

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
}

vector3<> randomUnitStep()
{
	vector3<> step(Random::uniform(), Random::uniform(), Random::uniform());
	step = normalize(step);
	return step;
}

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

vector3<> wrapPBC(vector3<> pos)
{
	for (int i = 0; i < 3; ++i)
	{
		if 		(pos[i] < 0) 		  pos[i] += boxSize[i];
		else if (pos[i] > boxSize[i]) pos[i] -= boxSize[i];
	}
}

vector3<> randomConePos(const std::vector<unitedAtom> polymerChains, const int a, const int b, const double distance)
{
	double theta, phi;
	vector3<> localX, localY, localZ, newPos;

	phi = (180 - 109.5)*M_PI/180;

	// calculate local coordinate system for the cone
	localZ = normalize(polymerChains[a].pos - polymerChains[b].pos);

	if (abs(localZ[0]) > 0) localX = vector3<>(1,0,0) - localZ * localZ.x();
	else					localX = vector3<>(0,1,0) - localZ * localZ.y();
	localX = normalize(localX);
	localY = normalize(cross(localZ, localX));

	theta = Random::uniform(0, 359)*M_PI/180;

	newPos = distance*sin(phi)*cos(theta)*localX
	       + distance*sin(phi)*sin(theta)*localY
		   + distance*cos(phi)		     *localZ;
	newPos += polymerChains[a].pos;
	return newPos;
}


bool checkCollision(
	const std::vector<unitedAtom> polymerChains,
	const std::vector<unitedAtom> chainLinks,
	const bool initFlag)
{
	bool flag = false;
	int iMax, newChainSize, chainSize;

	chainSize = polymerChains.size();
	newChainSize = chainLinks.size();

	iMax = chainSize - 1;
	if (initFlag) iMax = chainSize;
	// switch(newChainSize)
	// {
	// 	case 2: iMax = chainSize; 	  break;	// for CH3- CH at the initiation step
	// 	case 1:  break;	// for CH attached to the side branch
	// 	case 6:	iMax = chainSize - 1; break; 	// for just the side branch
	// 	case 7: iMax = chainSize - 1; break;	// for side branch with the next CH2
	// }
	
	for (int i = 0; i < iMax; ++i)
	{
		for (int j = 0; j < newChainSize; ++j)
		{
			double dist = ((polymerChains[i].pos) - (chainLinks[j].pos)).length();
			if (dist < minDist)
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
	FILE* fp = fopen("polymer.xyz", "w");
	fprintf(fp, "%d\n", polymerChains.size());
	fprintf(fp, "POLYSTYRENE\n");
	for (unitedAtom atom : polymerChains)
		fprintf(fp, "C\t%lf\t%lf\t%lf\n", atom.type, atom.pos[0], atom.pos[1], atom.pos[2]);
	fclose(fp);
}

void exportXSF(const std::vector<unitedAtom> polymerChains)
{
	FILE* fp = fopen("polymer.xsf", "w");
	fprintf(fp, "CRYSTAL\nPRIMVEC\n");
	fprintf(fp, "%lf\t0.0\t0.0\n0.0\t%lf\t0.0\n0.0\t0.0\t%lf\n", boxSize[0], boxSize[1], boxSize[2]);

	fprintf(fp, "CONVVEC\n");
	fprintf(fp, "%lf\t0.0\t0.0\n0.0\t%lf\t0.0\n0.0\t0.0\t%lf\n", boxSize[0], boxSize[1], boxSize[2]);

	fprintf(fp, "PRIMCOORD\n%d\t%d\n", polymerChains.size(), 1);
	for (unitedAtom atom : polymerChains)
		fprintf(fp, "C\t%lf\t%lf\t%lf\n", atom.pos[0], atom.pos[1], atom.pos[2]);
	fclose(fp);
}

void exportChainLengths()
{
	FILE* fp = fopen("chains.data", "w");
	for(int i=chainLengths.size()-1; i>0; i--) chainLengths[i] -= chainLengths[i-1];
	for(int i : chainLengths) fprintf(fp, "%d\n", i/8);
	fclose(fp);
}

// ring and CH2 together
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
	while(iTrial++ < maxTrials)
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

		if(checkCollision(polymerChains, newChainLinks)) continue;
		polymerChains.insert(polymerChains.end(), newChainLinks.begin(), newChainLinks.end());
		break;
	}

	if (iTrial < maxTrials) return true;
	
	iTrial = 0;
	while(iTrial++ < maxTrials)
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

		if(checkCollision(polymerChains, newChainLinks)) continue;
		polymerChains.insert(polymerChains.end(), newChainLinks.begin(), newChainLinks.end());
		break;
	}
	// printf("Tried %d times for just the ring at %d\n", iTrial, index);
	// if (iTrial >= maxTrials) return false;
	
	// iTrial = 0;
	// newChainLinks.resize(1);
	// while(iTrial++ < maxTrials)
	// {
	// 	ringPos = randomConePos(polymerChains, index-1, index-2, bondLengths[1]);
	// 	r12 = normalize(polymerChains[index-2].pos - rCH);
	// 	r10 = normalize(ringPos - rCH);

	// 	tetrahedral(r12, r10, pos1, pos2);
	// 	// discrete choice
	// 	rCH2 = (Random::uniform() > 0.5)? pos1 : pos2;
	// 	newChainLinks[0] = unitedAtom(2, bondLengths[1]*rCH2 + rCH);

	// 	if(checkCollision(polymerChains, newChainLinks)) continue;
	// 	polymerChains.insert(polymerChains.end(), newChainLinks.begin(), newChainLinks.end());
	// 	break;
	// }
	// printf("Tried %d times for CH2 after the ring at %d\n", iTrial, index);
	// if (iTrial < maxTrials) return true;
	return false;

	// add bonds data
	// listBonds.push_back(Bond(3, index, index+1));
	// for (int i = 0; i < 5; ++i)
	// 	listBonds.push_back(Bond(2, index+i+1, index+i+2));
	// listBonds.push_back(Bond(2, index+6, index+1));

	// listAngles.push_back(Angle(1, next_nearest_CH+1, index, index+1));
	// listAngles.push_back(Angle(2, index, index+1, index+2));
	// for (int i = 1; i <= 4; ++i)
	// 	listAngles.push_back(Angle(3, index+i, index+i+1, index+i+2));
	// listAngles.push_back(Angle(3, index+5, index+6, index+1));
	// listAngles.push_back(Angle(3, index+6, index+1, index+2));

}

// CH
bool addFirstHalf(std::vector<unitedAtom> &polymerChains)
{
	int chainSize, iTrial;

	chainSize = polymerChains.size();
	std::vector<unitedAtom> newCH(1, unitedAtom());
	newCH[0].type = 3;
	
	iTrial = 0;
	while(iTrial++ < maxTrials)
	{
		newCH[0].pos = randomConePos(polymerChains, chainSize-1, chainSize-8, bondLengths[0]);
		if(checkCollision(polymerChains, newCH)) continue;
		polymerChains.insert(polymerChains.end(), newCH.begin(), newCH.end());
		break;
	}
	if (iTrial >= maxTrials) return false;
	return true;
}

/*
Try adding CH3 and CH to initiate a new chain
if successful - return the index of the next position 
if failed - return 0, which means no more chains can be added
*/
bool initiateChain(std::vector<unitedAtom> &polymerChains)
{
	vector3<> step, pos0;
	std::vector<unitedAtom> newChainLinks(2, unitedAtom());

	int iTrial = 0;
	while(iTrial++ < maxTrials)
	{
		pos0 = vector3<>(
			Random::uniform(0, boxSize[0]),
			Random::uniform(0, boxSize[1]),
			Random::uniform(0, boxSize[2]));

		step = bondLengths[0] * randomUnitStep();

		newChainLinks[0] = unitedAtom(1, pos0);
		newChainLinks[1] = unitedAtom(3, pos0 + step);

		if (checkCollision(polymerChains, newChainLinks, true)) continue;
		polymerChains.push_back(newChainLinks[0]);
		polymerChains.push_back(newChainLinks[1]);
		break;
	}
	if (iTrial >= maxTrials) return false;

	// listBonds.push_back(Bond(1, index+1, index+2));
	return true;
}

void terminateChain(std::vector<unitedAtom> &polymerChains)
{
	int chainSize = polymerChains.size();
	int nChains = 1;

	printf("Terminated!! at %d\n", chainSize);
	chainLengths.push_back(chainSize);
	switch(polymerChains[chainSize-1].type)
	{
		case 1:	// can't add CH
		case 2:
			polymerChains.pop_back();
			break;
		case 5: // can't add the CH2
			polymerChains[chainSize-1].type = 2;
			break;
		case 3:	// can't add the ring and/or CH2
			polymerChains.pop_back();
			polymerChains.pop_back();
			if (chainSize > 2) polymerChains[chainSize-3].type = 2;
			break;
	}
}