#include "Random.h"
#include "matrix3.h"
#include "scalar.h"
#include "Euler.h"
#include "LatticeUtils.h"
#include <iostream>
#include <fstream>

// TODO
// 1. periodic boundaries during collision check, using periodic-lookup routine
// 2. book-keeping all bonds, etc

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

const int nMonomers = 1;	// number of monomers in the simulation box, calculate using desired density
// int nChains = 1;
// std::vector<int> chainLengths;
// int nAtoms = 7*nMonomers + nChains;

// std::vector<Bond> listBonds;
// std::vector<Angle> listAngles;
// std::vector<Dihedral> listDihedrals(nAtoms, Dihedral());
// std::vector<Dihedral> listImpropers(nAtoms, Dihedral());

// int nDihedrals=0, nImpropers=0;
const double minDist = 2;
const int maxTrials = 100;
const vector3<> boxSize(20, 20, 20);

// helper functions
vector3<> randomUnitStep();
void tetrahedral(const vector3<>, const vector3<>, vector3<>&, vector3<>&);
void wrapPBC(vector3<>&, const vector3<>);
vector3<> randomConePos(const std::vector<unitedAtom>, const int, const int, const double);

bool checkCollision(const std::vector<unitedAtom>, const std::vector<unitedAtom>);

// export functions
void exportXYZ(const std::vector<unitedAtom>);
// void exportXSF(const std::vector<unitedAtom>&);
// void exportLAMMPS(const std::vector<unitedAtom>&);
// void exportChainLengths();

// chain operations
bool initiateChain(std::vector<unitedAtom>&);
bool addSecondHalf(std::vector<unitedAtom>&, const std::vector<unitedAtom>);
bool addFirstHalf(std::vector<unitedAtom>&);
void terminateChain(std::vector<unitedAtom>&);

void initialise(std::vector<unitedAtom> &);

int main(int argc, char *argv[])
{
	std::vector<unitedAtom> polymerChains;
	std::vector<unitedAtom> aromaticRing;

	printf("===Self Avoiding Random Walk===\n");
	initialise(aromaticRing);
	
	int iMonomer = 0;

	if(initiateChain(polymerChains))
		while(iMonomer++ < nMonomers)
		{
			if(addSecondHalf(polymerChains, aromaticRing))
				if(addFirstHalf(polymerChains)) continue;
			terminateChain(polymerChains);

			if(initiateChain(polymerChains)) continue;
			break;
		}

	// export to XYZ file
	exportXYZ(polymerChains);
	// exportXSF(polymerChains);
	// export to LAMMPS data file
	// exportLAMMPS(polymerChains);

	// exportChainLengths(polymerChains);
	// printf("nChains : %d\n", nChains);

	return 0;
}

void initialise(std::vector<unitedAtom> &aromaticRing)
{
	// ========= Aromatic Ring =============
	double d  = bondLengths[2];
	double dz = d * cos(M_PI/3);
	double dx = d * sin(M_PI/3);

	aromaticRing.push_back( unitedAtom(4, vector3<>(-dx, 0,   dz)) );
	aromaticRing.push_back( unitedAtom(4, vector3<>(-dx, 0,   dz + d)) );
	aromaticRing.push_back( unitedAtom(4, vector3<>(  0, 0, 2*dz + d)) );
	aromaticRing.push_back( unitedAtom(4, vector3<>( dx, 0,   dz + d)) );
	aromaticRing.push_back( unitedAtom(4, vector3<>( dx, 0,   dz)) );
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

void wrapPBC(vector3<> &pos, const vector3<> boxSize)
{
	for (int i = 0; i < 3; ++i)
	{
		if 		(pos[i] < 0) pos[i] += boxSize[i];
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


bool checkCollision(const std::vector<unitedAtom> polymerChains, const std::vector<unitedAtom> chainLinks)
{
	// bool flag = false;
	// int iMax;
	
	// switch(num)
	// {
	// 	case 1: iMax = index - 7; break;
	// 	case 2: iMax = index - 1; break;
	// 	case 6: iMax = index - 2; break;
	// }
	
	// for (int i = 0; i < iMax; ++i)
	// {
	// 	for (int j = index; j < index+num; ++j)
	// 	{
	// 		double dist = (wrapPBC(polymerChains[i].pos) - wrapPBC(polymerChains[j].pos)).length();
	// 		if (dist < minDist)
	// 		{
	// 			flag = true;
	// 			break;
	// 		}
	// 	}
	// 	if (flag) break;
	// }
	return false;
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

// void exportXSF(const std::vector<unitedAtom> &polymerChains)
// {
// 	FILE* fp = fopen("polymer.xsf", "w");
// 	fprintf(fp, "CRYSTAL\nPRIMVEC\n");
// 	fprintf(fp, "%lf\t0.0\t0.0\n0.0\t%lf\t0.0\n0.0\t0.0\t%lf\n", boxSize[0], boxSize[1], boxSize[2]);

// 	fprintf(fp, "CONVVEC\n");
// 	fprintf(fp, "%lf\t0.0\t0.0\n0.0\t%lf\t0.0\n0.0\t0.0\t%lf\n", boxSize[0], boxSize[1], boxSize[2]);

// 	fprintf(fp, "PRIMCOORD\n%d\t%d\n", nAtoms, 1);
// 	for (unitedAtom atom : polymerChains)
// 		fprintf(fp, "6\t%lf\t%lf\t%lf\n", atom.pos[0], atom.pos[1], atom.pos[2]);
// 	fclose(fp);
// }

// void exportLAMMPS(const std::vector<unitedAtom> &polymerChains)
// {
// 	int nBonds = listBonds.size();
// 	int nAngles = listAngles.size();

// 	FILE* fp = fopen("polymer.data", "w");
// 	fprintf(fp, "POLYSTYRENE CHAINS\n\n");

// 	fprintf(fp, "%d\tatoms\n", nAtoms);
// 	fprintf(fp, "%d\tbonds\n", nBonds);
// 	fprintf(fp, "%d\tangles\n", nAngles);
// 	fprintf(fp, "%d\tdihedrals\n", nDihedrals);
// 	fprintf(fp, "%d\timpropers\n", nImpropers);
// 	fprintf(fp, "\n");
// 	fprintf(fp, "5\tatom types\n");
// 	fprintf(fp, "3\tbond types\n");
// 	fprintf(fp, "5\tangle types\n");
// 	fprintf(fp, "1\tdihedral types\n");
// 	fprintf(fp, "1\timproper types\n");
// 	fprintf(fp, "\n");
// 	fprintf(fp, "%lf\t%lf\txlo\txhi\n", -0.5, 0.5);
// 	fprintf(fp, "%lf\t%lf\tylo\tyhi\n", -0.5, 0.5);
// 	fprintf(fp, "%lf\t%lf\tzlo\tzhi\n", -0.5, 0.5);
// 	fprintf(fp, "\n");
// 	fprintf(fp, "Masses\n");
// 	for (int i=0;i<5;i++) fprintf(fp, "%d\t%lf\n", i+1, atomMass[i]);
// 	fprintf(fp, "\n");
// 	// ========= Coefficients ==============

// 	// ========= Atoms ==============
// 	fprintf(fp, "Atoms\n");
// 	for (int i=0; i<nAtoms; i++)
// 		fprintf(fp, "%d\t%d\t%lf\t%lf\t%lf\n", i+1, polymerChains[i].type, polymerChains[i].pos[0], polymerChains[i].pos[1], polymerChains[i].pos[2]);
// 	fprintf(fp, "\n");

// 	fprintf(fp, "Bonds\n");
// 	for (int i = 0; i < nBonds; i++)
// 	{
// 		fprintf(fp, "%d\t%d", i+1, listBonds[i].type);
// 		for (int j = 0; j < 2; ++j)	fprintf(fp, "\t%d", listBonds[i].species[j]);
// 		fprintf(fp, "\n");
// 	}
// 	fprintf(fp, "\n");
	
// 	fprintf(fp, "Angles\n");
// 	for (int i = 0; i < nAngles; i++)
// 	{
// 		fprintf(fp, "%d\t%d", i+1, listAngles[i].type);
// 		for (int j = 0; j < 3; ++j)	fprintf(fp, "\t%d", listAngles[i].species[j]);
// 		fprintf(fp, "\n");
// 	}
// 	fprintf(fp, "\n");
	
// 	fprintf(fp, "Dihedrals\n");
// 	for (int i = 0; i < nDihedrals; i++)
// 	{
// 		fprintf(fp, "%d\t%d\t", i+1, listDihedrals[i].type);
// 		for (int j = 0; j < 4; ++j)	fprintf(fp, "%d\t", listDihedrals[i].species[j]);
// 		fprintf(fp, "\n");
// 	}
// 	fprintf(fp, "\n");
	
// 	fprintf(fp, "Impropers\n");
// 	for (int i = 0; i < nImpropers; i++)
// 	{
// 		fprintf(fp, "%d\t%d\t", i+1, listImpropers[i].type);
// 		for (int j = 0; j < 4; ++j)	fprintf(fp, "%d\t", listImpropers[i].species[j]);
// 		fprintf(fp, "\n");
// 	}
// 	fprintf(fp, "\n");
// 	fclose(fp);
// }

// void exportChainLengths()
// {
// 	FILE* fp = fopen("chains.data", "w");
// 	chainLengths.push_back(nAtoms-1);
// 	for(int i=nChains-1; i>0; i--) chainLengths[i] -= chainLengths[i-1] + 1;
// 	for(int i : chainLengths) fprintf(fp, "%d\n", i/7);
// 	fclose(fp);
// }

// ring and CH2 together
bool addSecondHalf(std::vector<unitedAtom> &polymerChains, const std::vector<unitedAtom> aromaticRing)
{
	// initialize local variables
	vector3<> ringPos, r12, r10, rCH, rCH2, pos1, pos2;
	matrix3<> localCoords;
	int iTrial, index;

	index = polymerChains.size();
	rCH  = polymerChains[index-1].pos;
	
	// try adding aromatic ring and CH2
	std::vector<unitedAtom> newChainLinks(7, unitedAtom());
	iTrial = 0;
	while(iTrial++ < maxTrials)
	{
		ringPos = randomConePos(polymerChains, index-1, index-2, bondLengths[1]);
		newChainLinks[0] = unitedAtom(5, ringPos);

		// calculate the plane to place the ring 
		r12 = normalize(polymerChains[index-2].pos - rCH);
		r10 = normalize(ringPos - rCH);
		localCoords.set_row(2, r10);
		localCoords.set_row(0, normalize(r12 - r10 * dot(r10, r12)));
		localCoords.set_row(1, normalize(cross(localCoords.row(2), localCoords.row(0))));

		// add the remaining part of the ring
		for (int i = 1; i <= 5; ++i)
			newChainLinks[i] = unitedAtom(aromaticRing[i].type,	ringPos + aromaticRing[i].pos * localCoords);

		tetrahedral(r12, r10, pos1, pos2);
		// discrete choice
		rCH2 = (Random::uniform() > 0.5)? pos1 : pos2;
		newChainLinks[6] = unitedAtom(2, bondLengths[1]*rCH2 + rCH);

		if(checkCollision(polymerChains, newChainLinks)) continue;
		polymerChains.insert(polymerChains.end(), newChainLinks.begin(), newChainLinks.end());
		break;
	}

	// if (iTrial) printf("Tried %d times for ring at %d\n", iTrial, index);
	if (iTrial == maxTrials) return false;

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

	return true;
}

// CH
bool addFirstHalf(std::vector<unitedAtom> &polymerChains)
{
	int indexCH2, iTrial;

	indexCH2 = polymerChains.size()-1;
	std::vector<unitedAtom> newCH(1, unitedAtom());
	newCH[0].type = 3;
	
	iTrial = 0;
	while(iTrial++ < maxTrials)
	{
		newCH[0].pos = randomConePos(polymerChains, indexCH2, indexCH2-7, bondLengths[0]);
		if(checkCollision(polymerChains, newCH)) continue;
		polymerChains.insert(polymerChains.end(), newCH.begin(), newCH.end());
		break;
	}
	if (iTrial == maxTrials) return false;
	return true;
}


// 	do{
// 		prob = Random::uniform();

		
// 		if (!checkCollision(polymerChains, index, 1)) break;
// 		if (++iTrial == maxTrials) break;
// 	}while(true);

// 	// if (iTrial) printf("Tried %d times for CH at %d\n", iTrial, index);

// 	listBonds.push_back(Bond(1, index-6, index+1));
// 	listAngles.push_back(Angle(1, index-5, index-6, index+1));

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

		if (checkCollision(polymerChains, newChainLinks)) continue;
		polymerChains.push_back(newChainLinks[0]);
		polymerChains.push_back(newChainLinks[1]);
		break;
	}
	// if (iTrial) printf("Tried %d times to initiate\n", iTrial);
	if (iTrial == maxTrials) return false;

	// listBonds.push_back(Bond(1, index+1, index+2));
	return true;
}

void terminateChain(std::vector<unitedAtom> &polymerChains)
{
	int lastLink;
	lastLink = polymerChains.size() - 1;

	printf("Terminated!! at %d\n", lastLink);
	// switch(polymerChains[lastLink].type)
	// {
	// 	case 1:	// can't add CH
	// 	case 2:
	// 		break;
	// 	case 5: // can't add the CH2
	// 		for (int i = 0; i < 7; ++i)
	// 			polymerChains[lastLink-1+i] = unitedAtom();

	// 		if (polymerChains[lastLink-2].type == 1)
	// 		{
	// 			polymerChains[lastLink-2] = unitedAtom();
	// 			indexRestart = (index-2);
	// 		}
	// 		else
	// 		{
	// 			polymerChains[lastLink-8].type = 2;
	// 			indexRestart = (index-1);
	// 		}
	// 		break;
	// 	case 3:	// can't add the ring + CH2
	// 		polymerChains[lastLink] = unitedAtom();
	// 		if (polymerChains[lastLink-1].type == 1)
	// 		{
	// 			polymerChains[lastLink-1] = unitedAtom();
	// 			indexRestart = (index-1);
	// 		}
	// 		else
	// 		{
	// 			polymerChains[lastLink-7].type = 2;
	// 			indexRestart = index;
	// 		}
	// 		break;
	// }

	// chainLengths.push_back(indexRestart - 1);
	// nChains++; nAtoms++;
	// polymerChains.resize(nAtoms);
	// return indexRestart;
}