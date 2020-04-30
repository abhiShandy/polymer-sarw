#ifndef SARW_H
#define SARW_H

#define Avogadro 6.022e23
#define BondLength 1.54
#define AtomMass 14.0

#include <core/vector3.h>
#include <core/Util.h>
#include <PeriodicLookup.h>

// GLOBAL STRUCTS
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

class SARW
{
public:
    // #### Member Variables ####
    int nChains;
    vector3<> boxSize;
    double minDist;
    double targetMassDensity;
    int nGraftChains;
    string graftLoc;
    bool roundRobin;
    string boundary;
    int maxTrials;
    string polymer;
    bool logProgress, logSteps;

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
    int actualCount() const { return polymerChains.size(); }
    int vol() const { return boxSize[0]*boxSize[1]*boxSize[2]; }
    int targetCount() const { return nearbyint((vol() * targetMassDensity * Avogadro * 1e-24) / AtomMass); }
    double actualMassDensity()  const { return AtomMass * actualCount()/(Avogadro * vol() * 1e-24); }
    double targetNumberDensity()  const { return Avogadro * targetMassDensity / AtomMass; }
    double actualNumberDensity()  const { return Avogadro * actualMassDensity() / AtomMass; }

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
    {   polymerChains.push_back(atom);
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
