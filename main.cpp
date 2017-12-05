/*
*	Developed by Abhishek Shandilya
*	TODO
	1. periodic box
	2. book-keeping bonds
	3. adding CH2 after the side-branch
	4. break 2nd half into smaller modules
	5. use parameter file
	6. attain practical polymer density
*/


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
		0.10	//CH_aro - CH_aro
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

// std::vector<Bond> listBonds;
// std::vector<Angle> listAngles;
// std::vector<Dihedral> listDihedrals(nAtoms, Dihedral());
// std::vector<Dihedral> listImpropers(nAtoms, Dihedral());
int boxSize = 10;

// helper functions
vector3<> randomUnitStep();
void tetrahedral(const vector3<>, const vector3<>, vector3<>&, vector3<>&);
vector3<> wrapPBC(const vector3<>);
vector3<> randomConePos(const std::vector<unitedAtom>, const int, const int, const double);

bool checkCollision(const std::vector<unitedAtom>, const std::vector<unitedAtom>, const bool flag = false);

// export functions
void exportXYZ(const std::vector<unitedAtom>);
void exportXSF(const std::vector<unitedAtom>);
void exportLAMMPS(const std::vector<unitedAtom>);
void exportChainLengths(const std::vector<unitedAtom>);

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
	exportLAMMPS(polymerChains);

	exportChainLengths(polymerChains);

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

	boxSize = nearbyint(pow(nMonomers*(sideBranch.size()+2)/density, 0.3333));
	printf("BoxSize = %d\n", boxSize);
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
			vector3<> dist = (polymerChains[i].pos - chainLinks[j].pos);
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
	FILE* fp = fopen("polymer.xyz", "w");
	printf("nAtoms = %d\n", polymerChains.size());
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
	// int nBonds = listBonds.size();
	// int nAngles = listAngles.size();

	FILE* fp = fopen("polymer.data", "w");
	fprintf(fp, "POLYSTYRENE CHAINS\n\n");

	fprintf(fp, "%d\tatoms\n", polymerChains.size());
	// fprintf(fp, "%d\tbonds\n", nBonds);
	// fprintf(fp, "%d\tangles\n", nAngles);
	// fprintf(fp, "%d\tdihedrals\n", nDihedrals);
	// fprintf(fp, "%d\timpropers\n", nImpropers);
	fprintf(fp, "\n");
	fprintf(fp, "5\tatom types\n");
	fprintf(fp, "3\tbond types\n");
	fprintf(fp, "5\tangle types\n");
	fprintf(fp, "1\tdihedral types\n");
	fprintf(fp, "1\timproper types\n");
	fprintf(fp, "\n");
	fprintf(fp, "%lf\t%lf\txlo\txhi\n", -0.5, 0.5);
	fprintf(fp, "%lf\t%lf\tylo\tyhi\n", -0.5, 0.5);
	fprintf(fp, "%lf\t%lf\tzlo\tzhi\n", -0.5, 0.5);
	fprintf(fp, "\n");

	fprintf(fp, "Masses\n");
	for (int i=0;i<5;i++) fprintf(fp, "%d\t%lf\n", i+1, atomMass[i]);
	fprintf(fp, "\n");
	// ========= Coefficients ==============

	// ========= Atoms ==============
	fprintf(fp, "Atoms\n");
	int i = 0;
	for (unitedAtom ua : polymerChains)
		fprintf(fp, "%d\t%d\t%lf\t%lf\t%lf\n", ++i, ua.type, ua.pos[0], ua.pos[1], ua.pos[2]);
	fprintf(fp, "\n");

	fprintf(fp, "Bonds\n");
	// for (int i = 0; i < polymerChains.size()-1; i++)
	// {
		// currentType = polymerChains[i].type;
		// nextType = polymerChains[i+1].type;
		// fprintf(fp, "%d\t%d", i+1, listBonds[i].type);
		// for (int j = 0; j < 2; ++j)	fprintf(fp, "\t%d", listBonds[i].species[j]);
		// fprintf(fp, "\n");
	// }
	// i = 0;
	// for(Bond b : listBonds)
	// 	fprintf(fp, "%d\t%d\t%d\t%d\n", ++i, b.type, b.ua[0], b.ua[1]);
	// fprintf(fp, "\n");
	
	// fprintf(fp, "Angles\n");
	// for (int i = 0; i < nAngles; i++)
	// {
	// 	fprintf(fp, "%d\t%d", i+1, listAngles[i].type);
	// 	for (int j = 0; j < 3; ++j)	fprintf(fp, "\t%d", listAngles[i].species[j]);
	// 	fprintf(fp, "\n");
	// }
	// fprintf(fp, "\n");
	
	// fprintf(fp, "Dihedrals\n");
	// for (int i = 0; i < nDihedrals; i++)
	// {
	// 	fprintf(fp, "%d\t%d\t", i+1, listDihedrals[i].type);
	// 	for (int j = 0; j < 4; ++j)	fprintf(fp, "%d\t", listDihedrals[i].species[j]);
	// 	fprintf(fp, "\n");
	// }
	// fprintf(fp, "\n");
	
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

void exportChainLengths(const std::vector<unitedAtom> polymerChains)
{
	FILE* fp = fopen("chains.data", "w");
	int iAtom = 0, actual_nMonomers = 0;
	std::vector<int> chainLengths;
	
	for(unitedAtom i : polymerChains)
	{
		iAtom++;
		// this won't work for poly-propylene
		if (i.type == 1) chainLengths.push_back(iAtom);
		if (i.type == 5) actual_nMonomers++;
	}
	printf("nChains = %d\n", chainLengths.size());
	printf("Desired nMonomers = %d\n", nMonomers);
	printf("Actual nMonomers = %d\n", actual_nMonomers);
	printf("Desired number density = %.2f\n", density);
	printf("Actual mass density = %.2f\n", (104 * actual_nMonomers)/pow(boxSize, 3));

	for (int i = 0; i < chainLengths.size()-1; ++i)
		chainLengths[i] = (chainLengths[i+1] - chainLengths[i])/8;
	chainLengths[chainLengths.size()-1] = (polymerChains.size()+1 - chainLengths.back())/8;

	for(int i : chainLengths) fprintf(fp, "%d\n", i);
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
	if(polymerChains.size()==0)
	{	
		// vector3<> CH3_pos(0,0,0);
		vector3<> CH3_pos (0.5 * boxSize * vector3<>(1,1,1));
		vector3<> CH_pos ( CH3_pos + normalize(vector3<>(1,1,1)) * bondLengths[0] );
		polymerChains.push_back( unitedAtom(1, CH3_pos) );
		polymerChains.push_back( unitedAtom(3,  CH_pos));
	}

	else
	{

		int iTrial = 0;
		while(iTrial++ < maxTrials)
		{
			pos0 = vector3<>(
				Random::uniform(0, boxSize),
				Random::uniform(0, boxSize),
				Random::uniform(0, boxSize));


			step = bondLengths[0] * randomUnitStep();

			newChainLinks[0] = unitedAtom(1, pos0);
			newChainLinks[1] = unitedAtom(3, pos0 + step);

			if (checkCollision(polymerChains, newChainLinks, true)) continue;
			polymerChains.push_back(newChainLinks[0]);
			polymerChains.push_back(newChainLinks[1]);
			break;
		}
		if (iTrial >= maxTrials) return false;
	}

	// listBonds.push_back(Bond(1, index+1, index+2));
	return true;
}

void terminateChain(std::vector<unitedAtom> &polymerChains)
{
	int chainSize = polymerChains.size();

	// printf("Terminated!! at %d\n", chainSize);
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