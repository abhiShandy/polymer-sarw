#ifndef SARW_H
#define SARW_H

#define Avogadro 6.022e23

#include <core/vector3.h>
#include <core/Util.h>
#include <PeriodicLookup.h>

/////////////////////// GLOBAL VARIABLES & STRUCTS /////////////////////
const char species[6][10] = {"", "CH3", "CH2", "CH", "CH_aro", "C_aro"};

const double bondLengths[] =
{
	1.54,	//CH3    - CH2      or CH2 - CH
	1.51,	//CH     - C_aro
	1.40	//CH_aro - CH_aro
};

const double atomMass[] = {15, 14, 13, 13, 12};

struct unitedAtom
{
	int type; //1. CH3 2. CH2 3. CH 4. Si
	int chainID;
	vector3<> pos;

	unitedAtom() : type(0), chainID(0), pos(vector3<>(0,0,0)) {}
	unitedAtom(int type, vector3<> pos) : type(type), pos(pos) {}
	unitedAtom(int type, int chainID, vector3<> pos) : type(type), chainID(chainID), pos(pos) {}
};

struct Bond
{
	int type; // 1: CH3-CH2; 2: CH2-CH2
	int species[2];
	
	Bond(int t=0, int a=0, int b=0)
	{
		type = t;
		species[0] = a; species[1] = b;
	}
};

struct Angle
{
	int type; // 1. 109.5; 2. 120
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

//double mass_density = 0.92; // in g/cc
//double graftFraction = 0.20;
//int seed_radius = 10;
//double number_density = N_avogadro*mass_density/14.0; // units = per cc
// int sideLength = nearbyint(pow(10,8)*pow(nUnitedAtoms/number_density, 0.3333));
class SARW
{
public:
	// #### Member Variables ####
	string boundary;
	vector3<> boxSize;
	double graftFraction;
	bool logProgress, logSteps;
	int maxAtoms;
	int maxTrials;
	double minDist;
	int nChains;
	string polymer;
	bool roundRobin;
	double targetMassDensity;

	std::vector<Bond> listBonds;
	std::vector<Angle> listAngles;
	std::vector<Dihedral> listDihedrals;
	std::vector<vector3<>> listGrafts;
	std::vector<unitedAtom> polymerChains; // contains all finalised united atoms
	std::vector<int> lastIndices; // index of last united atom of each chain in polymerChains[]
	std::vector<int> penultimateIndices; // index of penulutimate united atom of each chain in polymerChains[]
	std::vector<int> chainLengths; // length of each polymer chain
	PeriodicLookup plook;

	// #### Constructor ####
	SARW(int nChains, vector3<> boxSize, double minDist)
	: nChains(nChains), boxSize(boxSize), minDist(minDist), plook(boxSize, minDist)
	{
		lastIndices.resize(nChains);
		penultimateIndices.resize(nChains);
		chainLengths.resize(nChains);
	}

	// #### Member Functions ####
	// - ## important quantities ##
	int nBonds() const {return (int)listBonds.size();}
	int nAngles() const {return (int)listAngles.size();}
	int nDihedrals() const {return (int)listDihedrals.size();}
	int nGrafts() const {return (int)listGrafts.size();}
	int targetCount() const { return nChains * maxAtoms; }
	int actualCount() const { return polymerChains.size(); }
	int vol() const { return boxSize[0]*boxSize[1]*boxSize[2]; }
	double actualMassDensity()  const { return ( 14.0*actualCount()+2.0*nChains)/(Avogadro * vol() * 1e-24); }
	double targetNumberDensity()  const { return Avogadro * targetMassDensity /14.0; }
	double actualNumberDensity()  const { return Avogadro * actualMassDensity() /14.0; }

	// - ## export functions ##
	void calcAnglesDihedrals();
	void exportXYZ() const;
	void exportXSF() const;
	void exportLAMMPS() const;
	void report() const;
	// - ## chain operations ##
	bool initiateChain();
	bool propagateChain();
	void terminateChain();
	bool randomSeed(std::vector<unitedAtom>& newChainLinks, const int iChain);
	bool checkCollision(const std::vector<unitedAtom>, const int ignoreIndex = -1);

	//Add atom to polymerChains and update plook accordingly
	inline void addAtom(const unitedAtom& atom)
	{	polymerChains.push_back(atom);
		plook.addPoint(atom.pos); //update plook
	}

	static InitParams initialize(int argc, char** argv, const char* description); //!< wrap initSystemCmdLine from JDFTx
	void setLogFlags(string logProgressFlag, string logStepsFlag);
	vector3<> randomUnitStep();
	vector3<> randomConePos(const int, const int, const double);
	void setMaxAtoms();
	void readGrafts(string fname);

};

#endif // SARW_H
