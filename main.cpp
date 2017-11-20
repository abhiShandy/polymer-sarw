#include "Random.h"
#include "matrix3.h"
#include "scalar.h"
#include "Euler.h"
#include <iostream>
#include <fstream>

// TODO
// 1. periodic boundaries
// 2. incomplete chains - use std-vector

const char species[6][10] = {"", "CH3", "CH2", "CH", "CH_aro", "C_aro"};
// const char species[6][10] = {"", "C", "H", "B", "F", "N"};
// const int unitedAtomTypes[] =
// 	{
// 		0,
// 		1, //CH3
// 		2, //CH2
// 		3, //CH
// 		4, //CH_aro
// 	};

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

struct unitedAtom
{
	int type;
	vector3<> pos;
	unitedAtom(int type, vector3<> pos) : type(type), pos(pos)
	{
	}
};

int nMonomers = 20;	// number of monomers in the simulation box, calculate using desired density
int nChains = 1;
int nAtoms = 7*nMonomers + nChains;
double minDist = 2;
int maxTrials = 100;
vector3<> boxSize(20, 20, 20); 

void printAtoms(unitedAtom*, int);
vector3<> randomStep();
bool checkCollision(unitedAtom*, int, int num = 1);
void exportXYZ(unitedAtom *, int);
void tetrahedral(vector3<>, vector3<>, vector3<>&, vector3<>&);
bool addFirstHalf(unitedAtom *, int);
bool addSecondHalf(unitedAtom *, int);
int initiate(unitedAtom *, int);
int terminate(unitedAtom *, int);

int main(int argc, char const *argv[])
{
	printf("Self Avoiding Random Walk\n");
	unitedAtom *listAtoms = (unitedAtom*) malloc(nAtoms * sizeof(unitedAtom));
	vector3<> step, newPos;
	int iAtom = 0;

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
		else iAtom += 6;

		// Step 4: add CH from last CH
		if (!addFirstHalf(listAtoms, iAtom))
		{
			// if can't find a spot, terminate the current chain
			iAtom = terminate(listAtoms, iAtom);
			iAtom = initiate(listAtoms, iAtom);
		}
		else iAtom++;
	}

	// ============ Termination ====================
	// find any incomplete chain, either remove them or complete them


	// export to XYZ file
	exportXYZ(listAtoms, nAtoms);
	// export to LAMMPS data file

	return 0;
}

void printAtoms(unitedAtom *listAtoms, int n)
{
	for (int i = 0; i < n; ++i)
	{
		printf("%d %d %.2f %.2f %.2f\n", i+1, listAtoms[i].type, listAtoms[i].pos[0], listAtoms[i].pos[1], listAtoms[i].pos[2]);
	}
}

vector3<> randomStep()
{
	vector3<> step;
	step[0] = Random::uniform(); step[1] = Random::uniform(); step[2] = Random::uniform();
	step = normalize(step);
	return step;
}

bool checkCollision(unitedAtom *listAtoms, int start, int num)
{
	bool flag = false;
	int iMax = start + num - 8;
	for (int i = 0; i < iMax; ++i)
	{
		for (int j = start; j < start+num; ++j)
		{
			double dist = (listAtoms[i].pos - listAtoms[j].pos).length();
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

void exportXYZ(unitedAtom *listAtoms, int nAtoms)
{
	FILE* fp = fopen("polymer.xyz", "w");
	fprintf(fp, "%d\n", nAtoms);
	fprintf(fp, "POLYSTYRENE\n");
	for (int i = 0; i < nAtoms; ++i)
		fprintf(fp, "C\t%lf\t%lf\t%lf\n", listAtoms[i].type, listAtoms[i].pos[0], listAtoms[i].pos[1], listAtoms[i].pos[2]);
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
bool addSecondHalf(unitedAtom *listAtoms, int start)
{
	double theta, phi, d, dx, dz;
	vector3<> localX, localY, localZ, newPos, r12, r10;
	matrix3<> localCoords;
	int iTrial = 0;
	
	int next_nearest_CH = start==2? start-2 : start-8;
	phi = 180-109.5;
	localZ = normalize(listAtoms[start-1].pos - listAtoms[next_nearest_CH].pos);
	if (abs(localZ[0]) > 0)
		localX = vector3<>(1,0,0) - localZ * localZ.x();
	else
		localX = vector3<>(0,1,0) - localZ * localZ.y();
	localX = normalize(localX);
	localY = normalize(cross(localZ, localX));
	
	unitedAtom *aromaticRing = (unitedAtom*)malloc(5 * sizeof(unitedAtom));

	d  = bondLengths[3];
	dz = d * cos(M_PI/3);
	dx = d * sin(M_PI/3);

	aromaticRing[0] = unitedAtom(4, vector3<>(-dx, 0,   dz));
	aromaticRing[1] = unitedAtom(4, vector3<>(-dx, 0,   dz + d));
	aromaticRing[2] = unitedAtom(4, vector3<>(  0, 0, 2*dz + d));
	aromaticRing[3] = unitedAtom(4, vector3<>( dx, 0,   dz + d));
	aromaticRing[4] = unitedAtom(4, vector3<>( dx, 0,   dz));

	do{
		theta = Random::uniform(0, 359);
		newPos = bondLengths[2]*sin(phi*M_PI/180)*cos(theta*M_PI/180)*localX
		       + bondLengths[2]*sin(phi*M_PI/180)*sin(theta*M_PI/180)*localY
			   + bondLengths[2]*cos(phi*M_PI/180)		   			 *localZ;
		newPos += listAtoms[start-1].pos;
		listAtoms[start] = unitedAtom(5, newPos);

		r12 = normalize(listAtoms[start].pos-listAtoms[start-1].pos);
		r10 = normalize(listAtoms[next_nearest_CH].pos-listAtoms[start-1].pos);
		localCoords.set_row(2, r12);
		localCoords.set_row(0, normalize(r10 - r12*dot(r10, r12)));
		localCoords.set_row(1, normalize(cross(localCoords.row(2), localCoords.row(0))));

		for (int i = 0; i < 5; ++i)
		{
			listAtoms[start+i+1].type = aromaticRing[i].type;
			listAtoms[start+i+1].pos  = listAtoms[start].pos + aromaticRing[i].pos * localCoords;
		}
		// PBC_wrap(listAtoms, start, 6);

		if (!checkCollision(listAtoms, start, 6)) break;
		if (++iTrial == maxTrials) break;
	}while(true);
	if (iTrial) printf("Tried %d times for ring at %d\n", iTrial, start);
	if (iTrial == maxTrials) return false;
	return true;
}

// CH
bool addFirstHalf(unitedAtom *listAtoms, int start)
{
	vector3<> r1, r12, r10, r8, r9;
	double prob;
	int iTrial = 0;
	
	r12 = listAtoms[start-6].pos-listAtoms[start-7].pos;
	r10 = listAtoms[start-8].pos-listAtoms[start-7].pos;
	r1  = listAtoms[start-7].pos;

	tetrahedral(r12, r10, r8, r9);

	do{
		prob = Random::uniform();
		r8 = (prob > 0.5)? r8 : r9;

		listAtoms[start] = unitedAtom(3, bondLengths[1]*r8+r1);
		// PBC_wrap(listAtoms, start, 1);
		
		if (!checkCollision(listAtoms, start, 1)) break;
		if (++iTrial == maxTrials) break;
	}while(true);

	if (iTrial) printf("Tried %d times for CH at %d\n", iTrial, start);
	if (iTrial == maxTrials) return false;
	return true;
}

int initiate(unitedAtom *listAtoms, int index)
{
	int iTrial = 0;
	vector3<> step, pos0;

	do{
		pos0 = vector3<>(
			Random::uniform(0, boxSize[0]),
			Random::uniform(0, boxSize[1]),
			Random::uniform(0, boxSize[2]));
		listAtoms[index] = unitedAtom(1, pos0);
		step  = bondLengths[0] * randomStep();
		step += pos0;
		listAtoms[index+1] = unitedAtom(3, step);
		// PBC_wrap(listAtoms, start, 2);

		if (!checkCollision(listAtoms, index, 2))	break;
		if (++iTrial == maxTrials) break;
	}while(true);
	if (iTrial) printf("Tried %d times to initiate at %d\n", iTrial, index);
	if (iTrial == maxTrials) return 0;
	else return (index+2);
}

int terminate(unitedAtom *listAtoms, int index)
{
	nChains++;
	printf("Terminated!! at %d\n", index);
	switch(listAtoms[index].type)
	{
		case 5:
			for (int i = 0; i < 7; ++i)
				listAtoms[index-1+i] = unitedAtom(0, vector3<>(0, 0, 0));

			if (listAtoms[index-2].type == 1)
			{
				listAtoms[index-2] = unitedAtom(0, vector3<>(0, 0, 0));
				return (index-2);
			}
			else
			{
				listAtoms[index-8].type = 2;
				return (index-1);
			}
			break;
		case 3:
			listAtoms[index] = unitedAtom(0, vector3<>(0, 0, 0));
			if (listAtoms[index-1].type == 1)
			{
				listAtoms[index-1] = unitedAtom(0, vector3<>(0, 0, 0));
				return (index-1);
			}
			else
			{
				listAtoms[index-7].type = 2;
				return index;
			}
			break;
	}
}

// // void PBC_wrap(unitedAtom *listAtoms, int start, int count)
// {
// 	// vector3<> pos;
// 	// for (int i = 0; i < count; ++i)
// 	// {
// 	// 	pos = listAtoms[start+i].pos;
// 	// 	if (pos[0] < 0) pos[0] += boxSize[0];
// 	// 	if (pos[1] < 0) pos[1] += boxSize[1];
// 	// 	if (pos[2] < 0) pos[2] += boxSize[2];
// 	// }
// }