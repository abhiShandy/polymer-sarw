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
		1.54,	//CH3 	 - CH2
		1.54,	//CH2 	 - CH
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

const int nMonomers = 10;	// number of monomers in the simulation box, calculate using desired density
int nChains = 1;
std::vector<int> chainLengths;
int nAtoms = 7*nMonomers + nChains;

std::vector<unitedAtom> listAtoms(nAtoms, unitedAtom());
std::vector<Bond> listBonds;
std::vector<Angle> listAngles;
std::vector<Dihedral> listDihedrals(nAtoms, Dihedral());
std::vector<Dihedral> listImpropers(nAtoms, Dihedral());

int nDihedrals=0, nImpropers=0;
const double minDist = 2;
const int maxTrials = 100;
const vector3<> boxSize(20, 20, 20);

vector3<> randomStep();
bool checkCollision(const std::vector<unitedAtom>&, const int, const int);
void exportXYZ(const std::vector<unitedAtom>&);
void exportXSF(const std::vector<unitedAtom>&);
void exportLAMMPS(const std::vector<unitedAtom>&);
void exportChainLengths();
void tetrahedral(vector3<>, vector3<>, vector3<>&, vector3<>&);
bool addFirstHalf(std::vector<unitedAtom>&, const int);
bool addSecondHalf(std::vector<unitedAtom>&, const int);
int initiate(std::vector<unitedAtom>&, const int);
int terminate(std::vector<unitedAtom>&, const int);
vector3<> wrapPBC(vector3<>);

int main(int argc, char *argv[])
{
	printf("Self Avoiding Random Walk\n");
	vector3<> step, newPos;
	int iAtom = 0, iMonomer = 0;

	// ============ Initiation ====================
	iAtom = initiate(listAtoms, iAtom);
	
	// ============ Propagation ====================
	while(iAtom < nAtoms && iAtom > 0)
	{
		// Step 3: add C_aro and construct aromatic ring, using a random angle
		if (!addSecondHalf(listAtoms, iAtom))
		{
			// if can't find a spot, undo Step 2 and 1, and terminate the current chain
			iAtom = terminate(listAtoms, iAtom);
			iAtom = initiate(listAtoms, iAtom);
			if (!iAtom) break;
			continue;
		}
		else
		{
			iAtom += 6;
			iMonomer++;
			if (iMonomer == nMonomers) break;
		}

		// Step 4: add CH from last CH
		if (!addFirstHalf(listAtoms, iAtom))
		{
			// if can't find a spot, terminate the current chain
			iAtom = terminate(listAtoms, iAtom);
			iAtom = initiate(listAtoms, iAtom);
		}
		else iAtom++;
	}

	// export to XYZ file
	exportXYZ(listAtoms);
	exportXSF(listAtoms);
	// export to LAMMPS data file
	exportLAMMPS(listAtoms);

	exportChainLengths();
	printf("nChains : %d\n", nChains);

	return 0;
}

vector3<> randomStep()
{
	vector3<> step(Random::uniform(), Random::uniform(), Random::uniform());
	step = normalize(step);
	return step;
}

bool checkCollision(std::vector<unitedAtom> &listAtoms, int index, int num)
{
	bool flag = false;
	int iMax;
	
	switch(num)
	{
		case 1: iMax = index - 7; break;
		case 2: iMax = index - 1; break;
		case 6: iMax = index - 2; break;
	}
	
	for (int i = 0; i < iMax; ++i)
	{
		for (int j = index; j < index+num; ++j)
		{
			double dist = (wrapPBC(listAtoms[i].pos) - wrapPBC(listAtoms[j].pos)).length();
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

vector3<> wrapPBC(vector3<> pos)
{
	for (int i = 0; i < 3; ++i)
	{
		if 		(pos[i] < 0) pos[i] += boxSize[i];
		else if (pos[i] > boxSize[i]) pos[i] -= boxSize[i];
	}
	return pos;
}

void exportXYZ(const std::vector<unitedAtom> &listAtoms)
{
	FILE* fp = fopen("polymer.xyz", "w");
	fprintf(fp, "%d\n", nAtoms);
	fprintf(fp, "POLYSTYRENE\n");
	for (unitedAtom atom : listAtoms)
		fprintf(fp, "C\t%lf\t%lf\t%lf\n", atom.type, atom.pos[0], atom.pos[1], atom.pos[2]);
	fclose(fp);
}

void exportXSF(const std::vector<unitedAtom> &listAtoms)
{
	FILE* fp = fopen("polymer.xsf", "w");
	fprintf(fp, "CRYSTAL\nPRIMVEC\n");
	fprintf(fp, "%lf\t0.0\t0.0\n0.0\t%lf\t0.0\n0.0\t0.0\t%lf\n", boxSize[0], boxSize[1], boxSize[2]);

	fprintf(fp, "CONVVEC\n");
	fprintf(fp, "%lf\t0.0\t0.0\n0.0\t%lf\t0.0\n0.0\t0.0\t%lf\n", boxSize[0], boxSize[1], boxSize[2]);

	fprintf(fp, "PRIMCOORD\n%d\t%d\n", nAtoms, 1);
	for (unitedAtom atom : listAtoms)
		fprintf(fp, "6\t%lf\t%lf\t%lf\n", atom.pos[0], atom.pos[1], atom.pos[2]);
	fclose(fp);
}

void exportLAMMPS(const std::vector<unitedAtom> &listAtoms)
{
	int nBonds = listBonds.size();
	int nAngles = listAngles.size();

	FILE* fp = fopen("polymer.data", "w");
	fprintf(fp, "POLYSTYRENE CHAINS\n\n");

	fprintf(fp, "%d\tatoms\n", nAtoms);
	fprintf(fp, "%d\tbonds\n", nBonds);
	fprintf(fp, "%d\tangles\n", nAngles);
	fprintf(fp, "%d\tdihedrals\n", nDihedrals);
	fprintf(fp, "%d\timpropers\n", nImpropers);
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
	for (int i=0; i<nAtoms; i++)
		fprintf(fp, "%d\t%d\t%lf\t%lf\t%lf\n", i+1, listAtoms[i].type, listAtoms[i].pos[0], listAtoms[i].pos[1], listAtoms[i].pos[2]);
	fprintf(fp, "\n");

	fprintf(fp, "Bonds\n");
	for (int i = 0; i < nBonds; i++)
	{
		fprintf(fp, "%d\t%d", i+1, listBonds[i].type);
		for (int j = 0; j < 2; ++j)	fprintf(fp, "\t%d", listBonds[i].species[j]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "Angles\n");
	for (int i = 0; i < nAngles; i++)
	{
		fprintf(fp, "%d\t%d", i+1, listAngles[i].type);
		for (int j = 0; j < 3; ++j)	fprintf(fp, "\t%d", listAngles[i].species[j]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "Dihedrals\n");
	for (int i = 0; i < nDihedrals; i++)
	{
		fprintf(fp, "%d\t%d\t", i+1, listDihedrals[i].type);
		for (int j = 0; j < 4; ++j)	fprintf(fp, "%d\t", listDihedrals[i].species[j]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "Impropers\n");
	for (int i = 0; i < nImpropers; i++)
	{
		fprintf(fp, "%d\t%d\t", i+1, listImpropers[i].type);
		for (int j = 0; j < 4; ++j)	fprintf(fp, "%d\t", listImpropers[i].species[j]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fclose(fp);
}

void exportChainLengths()
{
	FILE* fp = fopen("chains.data", "w");
	chainLengths.push_back(nAtoms-1);
	for(int i=nChains-1; i>0; i--) chainLengths[i] -= chainLengths[i-1] + 1;
	for(int i : chainLengths) fprintf(fp, "%d\n", i/7);
	fclose(fp);
}

void tetrahedral(vector3<> a, vector3<> b, vector3<> &c, vector3<> &d)
{
	a = normalize(a);
	b = normalize(b);

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

// C_aro and 5 CH_aro to form the aromatic ring
bool addSecondHalf(std::vector<unitedAtom> &listAtoms, int index)
{
	double theta, phi = 180 - 109.5;
	vector3<> localX, localY, localZ, newPos, r12, r10;
	matrix3<> localCoords;
	int iTrial = 0;

	// ========= Aromatic Ring =============
	std::vector<unitedAtom> aromaticRing(5, unitedAtom());

	double d  = bondLengths[3];
	double dz = d * cos(M_PI/3);
	double dx = d * sin(M_PI/3);

	aromaticRing[0] = unitedAtom(4, vector3<>(-dx, 0,   dz));
	aromaticRing[1] = unitedAtom(4, vector3<>(-dx, 0,   dz + d));
	aromaticRing[2] = unitedAtom(4, vector3<>(  0, 0, 2*dz + d));
	aromaticRing[3] = unitedAtom(4, vector3<>( dx, 0,   dz + d));
	aromaticRing[4] = unitedAtom(4, vector3<>( dx, 0,   dz));
	
	int next_nearest_CH = listAtoms[index-2].type == 1? index-2 : index-8;

	localZ = normalize(listAtoms[index-1].pos - listAtoms[next_nearest_CH].pos);

	if (abs(localZ[0]) > 0) localX = vector3<>(1,0,0) - localZ * localZ.x();
	else					localX = vector3<>(0,1,0) - localZ * localZ.y();
	localX = normalize(localX);
	localY = normalize(cross(localZ, localX));
	
	do{
		theta = Random::uniform(0, 359);
		newPos = bondLengths[2]*sin(phi*M_PI/180)*cos(theta*M_PI/180)*localX
		       + bondLengths[2]*sin(phi*M_PI/180)*sin(theta*M_PI/180)*localY
			   + bondLengths[2]*cos(phi*M_PI/180)		   			 *localZ;
		newPos += listAtoms[index-1].pos;
		listAtoms[index] = unitedAtom(5, newPos);

		r12 = normalize(listAtoms[index].pos-listAtoms[index-1].pos);
		r10 = normalize(listAtoms[next_nearest_CH].pos-listAtoms[index-1].pos);
		localCoords.set_row(2, r12);
		localCoords.set_row(0, normalize(r10 - r12*dot(r10, r12)));
		localCoords.set_row(1, normalize(cross(localCoords.row(2), localCoords.row(0))));

		for (int i = 0; i < 5; ++i)
		{
			listAtoms[index+i+1] = unitedAtom(aromaticRing[i].type,
				listAtoms[index].pos + aromaticRing[i].pos * localCoords);
		}


		if (!checkCollision(listAtoms, index, 6)) break;
		if (++iTrial == maxTrials) break;
	}while(true);

	// if (iTrial) printf("Tried %d times for ring at %d\n", iTrial, index);
	if (iTrial == maxTrials) return false;

	// add bonds data
	listBonds.push_back(Bond(3, index, index+1));
	for (int i = 0; i < 5; ++i)
		listBonds.push_back(Bond(2, index+i+1, index+i+2));
	listBonds.push_back(Bond(2, index+6, index+1));

	listAngles.push_back(Angle(1, next_nearest_CH+1, index, index+1));
	listAngles.push_back(Angle(2, index, index+1, index+2));
	for (int i = 1; i <= 4; ++i)
		listAngles.push_back(Angle(3, index+i, index+i+1, index+i+2));
	listAngles.push_back(Angle(3, index+5, index+6, index+1));
	listAngles.push_back(Angle(3, index+6, index+1, index+2));

	return true;
}

// CH
bool addFirstHalf(std::vector<unitedAtom> &listAtoms, int index)
{
	vector3<> r1, r12, r10, r8, r9;
	double prob;
	int iTrial = 0;
	
	r12 = listAtoms[index-6].pos-listAtoms[index-7].pos;
	r10 = listAtoms[index-8].pos-listAtoms[index-7].pos;
	r1  = listAtoms[index-7].pos;

	tetrahedral(r12, r10, r8, r9);

	do{
		prob = Random::uniform();
		r8 = (prob > 0.5)? r8 : r9;

		listAtoms[index] = unitedAtom(3, bondLengths[1]*r8+r1);
		
		if (!checkCollision(listAtoms, index, 1)) break;
		if (++iTrial == maxTrials) break;
	}while(true);

	// if (iTrial) printf("Tried %d times for CH at %d\n", iTrial, index);
	if (iTrial == maxTrials) return false;

	listBonds.push_back(Bond(1, index-6, index+1));
	listAngles.push_back(Angle(1, index-5, index-6, index+1));
	return true;
}

/*
Try adding CH3 and CH to initiate a new chain
if successful - return the index of the next position 
if failed - return 0, which means no more chains can be added
*/
int initiate(std::vector<unitedAtom> &listAtoms, int index)
{
	int iTrial = 0;
	vector3<> step, pos0;

	do{
		pos0 = vector3<>(
			Random::uniform(0, boxSize[0]),
			Random::uniform(0, boxSize[1]),
			Random::uniform(0, boxSize[2]));

		step = bondLengths[0] * randomStep();

		listAtoms[index]   = unitedAtom(1, pos0);
		listAtoms[index+1] = unitedAtom(3, pos0 + step);


		if (!checkCollision(listAtoms, index, 2))	break;
		if (++iTrial == maxTrials) break;
	}while(true);
	// if (iTrial) printf("Tried %d times to initiate at %d\n", iTrial, index);
	if (iTrial == maxTrials) return 0;

	listBonds.push_back(Bond(1, index+1, index+2));
	return (index+2);
}

int terminate(std::vector<unitedAtom> &listAtoms, int index)
{
	int indexRestart;

	printf("Terminated!! at %d\n", index);
	switch(listAtoms[index].type)
	{
		case 5:
			for (int i = 0; i < 7; ++i)
				listAtoms[index-1+i] = unitedAtom();

			if (listAtoms[index-2].type == 1)
			{
				listAtoms[index-2] = unitedAtom();
				indexRestart = (index-2);
			}
			else
			{
				listAtoms[index-8].type = 2;
				indexRestart = (index-1);
			}
			break;
		case 3:
			listAtoms[index] = unitedAtom();
			if (listAtoms[index-1].type == 1)
			{
				listAtoms[index-1] = unitedAtom();
				indexRestart = (index-1);
			}
			else
			{
				listAtoms[index-7].type = 2;
				indexRestart = index;
			}
			break;
	}

	chainLengths.push_back(indexRestart - 1);
	nChains++; nAtoms++;
	listAtoms.resize(nAtoms);
	return indexRestart;
}