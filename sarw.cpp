/*
*   Developed by Abhishek Shandilya
*   to-do
    1. periodic lookup
*/

// latest strategy - 29 Jan 2018 - construct the backbone and then add side-branches

#include <core/Random.h>
#include <core/matrix3.h>
#include <fluid/Euler.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "sarw.h"

/////////////////////// GLOBAL VARIABLES & STRUCTS /////////////////////
const char species[6][10] = {"", "CH3", "CH2", "CH", "CH_aro", "C_aro"};

const double bondLengths[] =
    {
        1.54,   //CH3    - CH2      or CH2 - CH
        1.51,   //CH     - C_aro
        1.40    //CH_aro - CH_aro
    };

const double atomMass[] = {15, 14, 13, 13, 12};

std::vector<vector3<>> grafts;

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

//////////////////////////////////// helper functions /////////////////////////////////////
vector3<> randomUnitStep();
vector3<> randomConePos(const std::vector<unitedAtom>, const int, const int, const double);

bool checkCollision(const std::vector<unitedAtom>, const std::vector<unitedAtom>, const int ignoreIndex = -1);

/////////////////////////////////// export functions //////////////////////////////////////
void exportXYZ(const std::vector<unitedAtom>);
void exportXSF(const std::vector<unitedAtom>);
void exportLAMMPS(const std::vector<unitedAtom>);
void printReport(const std::vector<unitedAtom>, std::vector<int>);

//////////////////////////////////// chain operations ////////////////////////////////////
bool initiateChain(std::vector<unitedAtom>&, std::vector<int>&, std::vector<int>&, std::vector<int> &);
bool propagateChain(std::vector<unitedAtom>&, std::vector<int>&, std::vector<int>&, std::vector<int> &);
void terminateChain(std::vector<unitedAtom>&, const std::vector<int>);
bool randomSeed(std::vector<unitedAtom>&, const std::vector<unitedAtom>, const int, const std::vector<vector3<>>);

int main(int argc, char *argv[])
{
    std::vector<unitedAtom> polymerChains;  // contains all finalised united atoms
    std::vector<int> lastIndices(nChains); // index of last united atom of each chain in polymerChains[]
    std::vector<int> penultimateIndices(nChains); // index of penulutimate united atom of each chain in polymerChains[]
    std::vector<int> chainLengths(nChains); // length of each polymer chain
    
    printf("===Self Avoiding Random Walk===\n");
    
    if (initiateChain(polymerChains, lastIndices, penultimateIndices, chainLengths)) // try initiating chains, if successful move on to next step
        while (polymerChains.size() < nUnitedAtoms) // if total number of united atoms is less than the target
            if (!propagateChain(polymerChains, lastIndices, penultimateIndices, chainLengths)) break; // propagate the chain, unless it fails
    terminateChain(polymerChains, lastIndices); // terminate the chains by changing the type of the last unitedAtom

    printf("\n===Exporting Files===\n");
    exportXYZ(polymerChains); // export to XYZ file for OVITO
    exportXSF(polymerChains); // export to XSF file for VESTA
    exportLAMMPS(polymerChains); // export to LAMMPS data file

    printReport(polymerChains, chainLengths);

    return 0;
}

vector3<> randomUnitStep()
{
    vector3<> step(Random::uniform(-1,1), Random::uniform(-1,1), Random::uniform(-1,1));
    step = normalize(step);
    return step;
}

vector3<> randomConePos(const std::vector<unitedAtom> polymerChains, const int lastIndex, const int penultimateIndex, const double distance)
{
    double theta, phi;
    vector3<> localX, localY, localZ, newPos;

    phi = (180 - 109.5)*M_PI/180; //tetrahedral angle

    // calculate local coordinate system for the cone
    localZ = normalize(polymerChains[lastIndex].pos - polymerChains[penultimateIndex].pos);

    if (abs(localZ[0]) > 0) localX = vector3<>(1,0,0) - localZ * localZ.x();
    else                    localX = vector3<>(0,1,0) - localZ * localZ.y();
    localX = normalize(localX);
    localY = normalize(cross(localZ, localX));

    theta = Random::uniform(0, 359)*M_PI/180;

    newPos = distance*sin(phi)*cos(theta)*localX
           + distance*sin(phi)*sin(theta)*localY
           + distance*cos(phi)           *localZ;
    newPos += polymerChains[lastIndex].pos;
    return newPos;
}

/*
Check Collision of a potential section of chain
- return true if colliding
- return false if the chain can be added
*/
bool checkCollision(
    const std::vector<unitedAtom> polymerChains,
    const std::vector<unitedAtom> newChainLinks,
    const int ignoreIndex)
{
    int newChainSize, chainSize;

    chainSize    = polymerChains.size();
    newChainSize = newChainLinks.size();

    double zMin = 0.;
    for (vector3<> g : grafts) zMin = zMin < g.z() ? g.z() : zMin;

    for (int i = 0; i < chainSize; ++i)
    {
        if (ignoreIndex >= 0 && i==ignoreIndex) continue;
        for (int j = 0; j < newChainSize; ++j)
        {
            vector3<> dist = (polymerChains[i].pos - newChainLinks[j].pos);
            for (int k = 0; k < 3; ++k) // minimum image convention
                dist[k] -= boxSize[k] * nearbyint(dist[k] / boxSize[k]);
            if ((dist.length() < minDist) || (newChainLinks[j].pos[0] < 0) || (newChainLinks[j].pos[1] < 0) || (newChainLinks[j].pos[2] < zMin))
                return true;
        }
    }
    return false;
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

/*
Export DAT file for LAMMPS input
- TODO : use Class/Struct for passing all data structures
*/
std::vector<Bond> listBonds;
void exportLAMMPS(const std::vector<unitedAtom> polymerChains)
{
    printf("Exporting LAMMPS data file: polymer.dat\n");
    std::vector<Angle>      listAngles;
    std::vector<Dihedral>   listDihedrals;

    // ------------ Counting Angles and Dihedrals ----------------
    int a1, a2, a3, a4;
    int i = 0;
    for (int iChain = 1; iChain <= nChains; ++iChain)
    {
        a1=-1, a2=-1, a3=-1, a4=-1;
        // the very first monomer of iChain
        i = 2*iChain - 2; a4 = i;
        while(i < polymerChains.size())
        {
            // if all three indices are found, add to angles list
            if (a2>-1) listAngles.push_back(Angle(1, a2+1, a3+1, a4+1));
            // if all four indices are found, add to dihedrals list
            if(a1>-1) listDihedrals.push_back(Dihedral(1, a1+1, a2+1, a3+1, a4+1));
            // search for next set, and udpate all indices
            while(++i < polymerChains.size())
                if (polymerChains[i].chainID == iChain)
                {
                    a1 = a2; a2 = a3; a3 = a4; a4=i;
                    break;
                }
        }
    }
    // ------------ Counting impropers ----------------

    FILE* fp = fopen("polymer.dat", "w");
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
        for (int j = 0; j < 4; ++j) fprintf(fp, "%d\t", listDihedrals[i].species[j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    
    // fprintf(fp, "Impropers\n");
    // for (int i = 0; i < nImpropers; i++)
    // {
    //  fprintf(fp, "%d\t%d\t", i+1, listImpropers[i].type);
    //  for (int j = 0; j < 4; ++j) fprintf(fp, "%d\t", listImpropers[i].species[j]);
    //  fprintf(fp, "\n");
    // }
    // fprintf(fp, "\n");
    fclose(fp);
}

void printReport(const std::vector<unitedAtom> polymerChains, const std::vector<int> chainLengths)
{
    printf("\n===Report===\n");
    
    printf("\nBoxSize\t: %d x %d x %d\n", boxSize[0], boxSize[1], boxSize[2]);
    
    printf("\nnUnitedAtoms::\n");
    printf("Desired\t: %d\n", nUnitedAtoms);
    printf("Actual\t: %d\n", polymerChains.size());
    
    int vol = (boxSize[0]*boxSize[1]*boxSize[2]);
    printf("\nNumber Density::\n");
    printf("Desired\t: %.3e\n", number_density);
    printf("Actual\t: %.3e\n", polymerChains.size()/(vol*pow(10,-24)));

    printf("\nMass Density::\n");
    printf("Desired\t: %.2f\n", mass_density);
    printf("Actual\t: %.2f\n", ( 14*polymerChains.size() )/(N_avogadro*vol*pow(10,-24)) );

    printf("Exporting distribution of chain lengths: chains.dat\n");
    FILE* fp = fopen("chains.dat", "w");
    for (int i : chainLengths) fprintf(fp, "%d\n", i);
    fclose(fp);
}

/* Chain propagation
 - add CH2 on a cone of tetrahedral angle
 - TODO: add side-chain(s)
 */
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
        if (ROUND_ROBIN) iChain = i;
        else iChain = nearbyint(Random::uniform(0, nChains-1));

        lastIndex = lastIndices[iChain];
        penultimateIndex = penultimateIndices[iChain];
        
        int iTrial = 0;
        while (iTrial++ < maxTrials)
        {
            // avoid hard-coding any number
            newCH2[0] = unitedAtom(2, iChain+1, randomConePos(polymerChains, lastIndex, penultimateIndex, bondLengths[0]));
            if (checkCollision(polymerChains, newCH2, lastIndex)) continue;
            
            polymerChains.push_back(newCH2[0]);

            penultimateIndices[iChain] = lastIndices[iChain];
            lastIndices[iChain]        = polymerChains.size() - 1;
            listBonds.push_back(Bond(1, penultimateIndices[iChain] + 1, lastIndices[iChain] + 1));
            chainLengths[iChain]++;

            break;
        }
        flag = flag && (iTrial-1 >= maxTrials); // TODO: recheck the logic

        if (DEBUG) printf("DEBUG:: iTrial %d\n", iTrial);
        if (iTrial-1 < maxTrials)
        {
            if (DEBUG) printf("DEBUG:: Propagated %dth chain at %d\n", i, polymerChains.size());
            if (LOG) printf("LOG:: Progress = %d\n", 100*polymerChains.size()/nUnitedAtoms);
        }
    }
    return !flag;
}

/*
Try adding a pair of CH3 and CH2 to initiate a new chain
    if successful - add the united-atoms to data-structure, and return true
    if failed - return false, which means no more chains can be added
Input
    - file containing location of silica graft points (FILE FORMAT: x y z)


*/
bool initiateChain(std::vector<unitedAtom> &polymerChains, std::vector<int> &lastIndices, std::vector<int> &penultimateIndices, std::vector<int> &chainLengths)
{
    // Store graft locations in matrix
    std::ifstream graft_pos("graft.dat");
    float x,y,z;
    while (graft_pos >> x >> y >> z) {
        grafts.push_back(vector3<>(x,y,z));
    }
    
    if (DEBUG) printf("DEBUG:: Entered Initiation step\n");

    std::vector<unitedAtom> newChainLinks(2, unitedAtom());

    // Loop through each chain and initiate
    for (int iChain = 0; iChain < nChains; ++iChain)
    {
        if (randomSeed(newChainLinks, polymerChains, iChain, grafts))
        {
            newChainLinks[0].chainID = iChain+1;
            newChainLinks[1].chainID = iChain+1;
            polymerChains.push_back(newChainLinks[0]);
            polymerChains.push_back(newChainLinks[1]);
            // update indices, index starts from 0
            lastIndices[iChain] = polymerChains.size() - 1;
            penultimateIndices[iChain] = polymerChains.size() - 2;
            // update listBond, index starts from 1
            listBonds.push_back(Bond(1, penultimateIndices[iChain]+1, lastIndices[iChain]+1));
            chainLengths[iChain] += 2;
            if (DEBUG) printf("DEBUG:: Initiated %dth seed\n", iChain);
        }
        else
            return false;
    }
    if (DEBUG) printf("DEBUG:: Initiated %d seeds\n", nChains);

    return true;
}

/*
Place randomly located seeds:
Inputs:
    - an empty vector of 2 unitedAtoms
    - list of all unitedAtoms added so far
Outputs:
    - boolean: true if it successfully found a random location for the seed, otherwise false
    - location of random seed in newChainLinks
*/
bool randomSeed(std::vector<unitedAtom> &newChainLinks, 
    const std::vector<unitedAtom> polymerChains,
    const int iChain,
    const std::vector<vector3<>> grafts)
{
    vector3<> step, pos0;
    int iTrial = 0;
    while (iTrial++ < maxTrials)
    {
        if (iChain < grafts.size()) // if atom is grafted, seed from here
            pos0 = vector3<>(grafts[iChain][0], grafts[iChain][1], grafts[iChain][2]);
        else // if atom is not grafted, seed randomly
        {
            pos0 = vector3<>(
                Random::uniform(0, boxSize[0]),
                Random::uniform(0, boxSize[1]),
                Random::uniform(0, boxSize[2]));
        }

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
        polymerChains[lastIndices[i]].type = 1;
        if (DEBUG) printf("Terminated chain # %d\n", i+1);
    }
}