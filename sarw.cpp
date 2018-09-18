/*
*   Developed by Abhishek Shandilya
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <core/Random.h>
#include "InputMap.h"
#include "sarw.h"
#include "config.h"

int main(int argc, char **argv)
{
    InitParams ip = SARW::initialize(argc, argv, "Self Avoiding Random Walk");
    InputMap inputMap(ip.inputFilename);

    // set parameters from input file
    int nChains = inputMap.get("nChains", 1);
    SARW s = SARW(nChains);
    s.minChainLength = inputMap.get("minChainLength", 2);
    s.boxSize = inputMap.getVector("boxSize", vector3<>(10,10,10));
    s.maxTrials = inputMap.get("maxTrials", 100);
    s.targetMassDensity = inputMap.get("targetMassDensity", 0.92);
    s.minDist = inputMap.get("minDist", 2.0);
    s.roundRobin = inputMap.getString("roundRobin")=="yes";
    s.polymer = inputMap.getString("polymer");
    s.setLogFlags(inputMap.getString("logProgress"), inputMap.getString("logSteps"));
    s.graftedSeeds = inputMap.getString("graftedSeeds");
    s.graftFraction = inputMap.get("graftFraction", 0.);
    s.growthBias = inputMap.getVector("growthBias", vector3<>(0,0,0));
    s.boundary = inputMap.getString("boundary");

    logPrintf("\nINPUTS:\n");
    logPrintf("nChains = %d\n", s.nChains);
    logPrintf("minChainLength = %d\n", s.minChainLength);
    logPrintf("maxTrials = %d\n", s.maxTrials);
    logPrintf("boxSize = (%lg x %lg x %lg)\n", s.boxSize[0], s.boxSize[1], s.boxSize[2]);
    logPrintf("targetMassDensity = %lg\n", s.targetMassDensity);
    logPrintf("minDist = %lg\n", s.minDist);
    logPrintf("polymer = %s\n", s.polymer.c_str());

    logPrintf("\nFLAGS:\n");
    logPrintf("roundRobin = %s\n", s.roundRobin?"true":"false");
    logPrintf("logProgress = %s\n", s.logProgress?"true":"false");
    logPrintf("logSteps = %s\n", s.logSteps?"true":"false");

    logPrintf("\nCONSTRAINTS:\n");
    logPrintf("graftedSeeds = %s\n", s.graftedSeeds.c_str());
    logPrintf("growthBias = (%lg x %lg x %lg)\n", s.growthBias[0], s.growthBias[1], s.growthBias[2]);
    logPrintf("boundary = %s\n", s.boundary.c_str());

    if (s.logProgress) logPrintf("\nPROGRESS: ");
    if (s.initiateChain()) // try initiating chains, if successful move on to next step
        while (s.actualCount() < s.targetCount()) // if total number of united atoms is less than the target
        {
            if (!s.propagateChain()) break; // propagate the chain, unless it fails
            if (s.logProgress) logPrintf("%d, ", 100*s.actualCount()/s.targetCount());
        }
    s.terminateChain(); // terminate the chains by changing the type of the last unitedAtom
    if (s.logProgress) logPrintf("\n");

    s.calcAnglesDihedrals();
    s.exportLAMMPS();
    s.exportXYZ();
    s.exportXSF();
    s.report();

    return 0;
}

InitParams SARW::initialize(int argc, char** argv, const char* description)
{   InitParams ip;
    ip.packageName = PACKAGE_NAME;
    ip.versionString = VERSION_STRING;
    ip.versionHash = GIT_HASH;
    ip.description = description;
    initSystemCmdline(argc, argv, ip);
    // fftw_mpi_init();
    return ip;
}

vector3<> SARW::randomUnitStep()
{
    vector3<> step(Random::uniform(-.5,.5), Random::uniform(-.5,.5), Random::uniform(-.5,.5));

    if (growthBias[0]==1) step[0] = Random::uniform(0,1);
    else if (growthBias[0]==-1) step[0] = Random::uniform(-1,0);

    if (growthBias[1]==1) step[1] = Random::uniform(0,1);
    else if (growthBias[1]==-1) step[1] = Random::uniform(-1,0);

    if (growthBias[2]==1) step[2] = Random::uniform(0,1);
    else if (growthBias[2]==-1) step[2] = Random::uniform(-1,0);
    
    return normalize(step);
}

/*
 * Generate a random point on a Cone of Tetrahedral head angle given the following inputs:
 * - lastIndex, penultimateIndex : indices of unitedAtoms to form the axis of the cone
 * - distance: distance along the cone's surface (bond length between the last unitedAtom and the newly created one)
 */
vector3<> SARW::randomConePos(const int lastIndex, const int penultimateIndex, const double distance)
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
  
    while(true)
    {
      theta = Random::uniform(0, 359)*M_PI/180;

      newPos = distance*sin(phi)*cos(theta)*localX
	     + distance*sin(phi)*sin(theta)*localY
	     + distance*cos(phi)           *localZ;
      if (growthBias[0] && newPos[0]>0) break;
      if (growthBias[1] && newPos[1]>0) break;
      if (growthBias[2] && newPos[2]>0) break;
      if (!(growthBias[0]+growthBias[1]+growthBias[2])) break;
    }
    
    newPos += polymerChains[lastIndex].pos;
    return newPos;
}

/*
Check Collision of a potential section of chain
- return true if colliding
- return false if the chain can be added
*/
bool SARW::checkCollision(const std::vector<unitedAtom> newChainLinks, const int ignoreIndex)
{
    int chainSize    = actualCount();
    int newChainSize = newChainLinks.size();

    // Check for collisions
    for (int i = 0; i < chainSize; ++i)
    {
        if (ignoreIndex >= 0 && i==ignoreIndex) continue;
        for (int j = 0; j < newChainSize; ++j)
        {
            vector3<> dist = (polymerChains[i].pos - newChainLinks[j].pos);
            for (int k = 0; k < 3; ++k) // minimum image convention
	      dist[k] -= boxSize[k] * nearbyint(dist[k] / boxSize[k]);

            bool collisionFlag = (dist.length() < minDist);
            if (collisionFlag) return true;
        }
    }
    
    // Fixed boundary condtions:
    for (int j = 0; j < newChainSize; ++j)
    {
	bool collisionFlag = ((boundary=="ppf") && ((newChainLinks[j].pos[2]<0) || (newChainLinks[j].pos[2]>boxSize[2])));
	collisionFlag = collisionFlag || ((boundary=="pfp") && ((newChainLinks[j].pos[1]<0) || (newChainLinks[j].pos[1]>boxSize[1])));
	collisionFlag = collisionFlag || ((boundary=="fpp") && ((newChainLinks[j].pos[0]<0) || (newChainLinks[j].pos[0]>boxSize[0])));
	if (collisionFlag) return true;
    }
    return false;
}

/* Chain propagation
 - add CH2 on a cone of tetrahedral angle
 - TODO: add side-chain(s)
 */
bool SARW::propagateChain()
{
    if (logSteps) logPrintf("Entered Propagation step\n");
    // propagate all chains, if possible
    // udpate indices and lengths
    int iChain, lastIndex, penultimateIndex;
    bool flag = true;

    std::vector<unitedAtom> newCH2(1, unitedAtom());
    
    for (int i = 0; i < nChains; ++i)
    {
        if (roundRobin) iChain = i;
        else iChain = nearbyint(Random::uniform(0, nChains-1));

        lastIndex = lastIndices[iChain];
        penultimateIndex = penultimateIndices[iChain];
        
        int iTrial = 0;
        while (iTrial++ < maxTrials)
        {
            // avoid hard-coding any number
            newCH2[0] = unitedAtom(2, iChain+1, randomConePos(lastIndex, penultimateIndex, bondLengths[0]));
            if (checkCollision(newCH2, lastIndex)) continue;
            
            polymerChains.push_back(newCH2[0]);

            penultimateIndices[iChain] = lastIndices[iChain];
            lastIndices[iChain]        = actualCount() - 1;
            listBonds.push_back(Bond(1, penultimateIndices[iChain] + 1, lastIndices[iChain] + 1));
            chainLengths[iChain]++;

            break;
        }
        flag = flag && (iTrial-1 >= maxTrials); // TODO: recheck the logic

        if (logSteps) logPrintf("iTrial %d\n", iTrial);
        if ((iTrial-1 < maxTrials) && logSteps) logPrintf("Propagated %dth chain at %d\n", i, actualCount());
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
bool SARW::initiateChain()
{
    // Store graft locations in matrix
    if (graftedSeeds=="file")
    {
	logPrintf("Reading graft.dat\n");
        std::ifstream graft_pos("graft.dat");
        double x,y,z;
        while (graft_pos >> x >> y >> z) {
            listGrafts.push_back(vector3<>(x,y,z));
        }
    }
    
    if (logSteps) printf("Entered Initiation step\n");

    std::vector<unitedAtom> newChainLinks(2, unitedAtom());

    // Loop through each chain and initiate
    for (int iChain = 0; iChain < nChains; ++iChain)
    {
        if (randomSeed(newChainLinks, iChain))
        {
            newChainLinks[0].chainID = iChain+1;
            newChainLinks[1].chainID = iChain+1;
            polymerChains.push_back(newChainLinks[0]);
            polymerChains.push_back(newChainLinks[1]);
            // update indices, index starts from 0
            lastIndices[iChain] = actualCount() - 1;
            penultimateIndices[iChain] = actualCount() - 2;
            // update listBond, index starts from 1
            listBonds.push_back(Bond(1, penultimateIndices[iChain]+1, lastIndices[iChain]+1));
            chainLengths[iChain] += 2;
            if (logSteps) logPrintf("Initiated %dth seed\n", iChain);
        }
        else
            return false;
    }
    if (logSteps) logPrintf("Initiated %d seeds\n", nChains);

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
bool SARW::randomSeed(std::vector<unitedAtom>& newChainLinks, const int iChain)
{
    vector3<> step, pos0;
    int iTrial = 0;
    while (iTrial++ < maxTrials)
    {
        if (iChain < nGrafts()) // if atom is grafted, seed from here
            pos0 = vector3<>(listGrafts[iChain][0], listGrafts[iChain][1], listGrafts[iChain][2]);
        else // if atom is not grafted, seed randomly
        {   pos0 = vector3<>(
                Random::uniform(0, boxSize[0]),
                Random::uniform(0, boxSize[1]),
                Random::uniform(0, boxSize[2]));
        }
        if (iChain < nChains*graftFraction)
        {   if (graftedSeeds == "x") pos0[0]=0.;
            if (graftedSeeds == "y") pos0[1]=0.;
            if (graftedSeeds == "z") pos0[2]=0.;
        }

        // CH atom at a random position on a sphere around the CH3 seed
        step = bondLengths[0] * randomUnitStep();
        newChainLinks[0] = (unitedAtom(1, pos0));
        newChainLinks[1] = (unitedAtom(2, pos0 + step));
        if (!checkCollision(newChainLinks)) break;
    }
    if (iTrial >= maxTrials) return false;
    return true;
}

void SARW::terminateChain()
{
    if (logSteps) printf("Terminated!! at %d\n", actualCount());
    for (int i = 0; i < nChains; ++i)
    {
        polymerChains[lastIndices[i]].type = 1;
        if (logSteps) printf("Terminated chain # %d\n", i+1);
    }
}

void SARW::exportXYZ() const
{
    logPrintf("Exporting XYZ file: polymer.xyz\n");
    FILE* fp = fopen("polymer.xyz", "w");
    fprintf(fp, "%d\n", actualCount());
    fprintf(fp, "POLYSTYRENE\n");
    for (unitedAtom atom : polymerChains)
        fprintf(fp, "C\t%lf\t%lf\t%lf\n", atom.pos[0], atom.pos[1], atom.pos[2]);
    fclose(fp);
}

void SARW::exportXSF() const
{
    logPrintf("Exporting XSF file: polymer.xsf\n");
    FILE* fp = fopen("polymer.xsf", "w");
    fprintf(fp, "CRYSTAL\nPRIMVEC\n");
    fprintf(fp, "%lf\t0.0\t0.0\n0.0\t%lf\t0.0\n0.0\t0.0\t%lf\n", boxSize[0], boxSize[1], boxSize[2]);

    fprintf(fp, "CONVVEC\n");
    fprintf(fp, "%lf\t0.0\t0.0\n0.0\t%lf\t0.0\n0.0\t0.0\t%lf\n", boxSize[0], boxSize[1], boxSize[2]);

    fprintf(fp, "PRIMCOORD\n%d\t%d\n", actualCount(), 1);
    for (unitedAtom atom : polymerChains)
        fprintf(fp, "C\t%lf\t%lf\t%lf\n", atom.pos[0], atom.pos[1], atom.pos[2]);
    fclose(fp);
}

void SARW::calcAnglesDihedrals()
{
    // ------------ Counting Angles and Dihedrals ----------------
    int a1, a2, a3, a4;
    int i = 0;
    for (int iChain = 1; iChain <= nChains; ++iChain)
    {
        a1=-1, a2=-1, a3=-1, a4=-1;
        // the very first monomer of iChain
        i = 2*iChain - 2; a4 = i;
        while(i < actualCount())
        {
            // if all three indices are found, add to angles list
            if (a2>-1) listAngles.push_back(Angle(1, a2+1, a3+1, a4+1));
            // if all four indices are found, add to dihedrals list
            if(a1>-1) listDihedrals.push_back(Dihedral(1, a1+1, a2+1, a3+1, a4+1));
            // search for next set, and udpate all indices
            while(++i < actualCount())
                if (polymerChains[i].chainID == iChain)
                {
                    a1 = a2; a2 = a3; a3 = a4; a4=i;
                    break;
                }
        }
    }
}

/*
Export DAT file for LAMMPS input
*/
void SARW::exportLAMMPS() const
{
    logPrintf("\nExporting LAMMPS data file: polymer.dat\n");

    // ------------ Counting impropers ----------------

    FILE* fp = fopen("polymer.dat", "w");
    fprintf(fp, "%s chains\n\n", polymer.c_str());

    fprintf(fp, "%5d\tatoms\n", actualCount());
    fprintf(fp, "%5d\tbonds\n", nBonds());
    fprintf(fp, "%5d\tangles\n", nAngles());
    fprintf(fp, "%5d\tdihedrals\n", nDihedrals());
    // fprintf(fp, "%d\timpropers\n", nImpropers);
    fprintf(fp, "\n");
    fprintf(fp, "2\tatom types\n");
    fprintf(fp, "1\tbond types\n");
    fprintf(fp, "1\tangle types\n");
    fprintf(fp, "1\tdihedral types\n");
    // fprintf(fp, "1\timproper types\n");
    fprintf(fp, "\n");
    fprintf(fp, "0.0 %lf xlo xhi\n", boxSize[0]);
    fprintf(fp, "0.0 %lf ylo yhi\n", boxSize[1]);
    fprintf(fp, "0.0 %lf zlo zhi\n", boxSize[2]);
    fprintf(fp, "\n");

    fprintf(fp, "Masses\n\n");
    for (int i=0;i<2;i++) fprintf(fp, "%d\t%lf\n", i+1, atomMass[i]);
    fprintf(fp, "\n");

    // ========= Atoms ==============
    fprintf(fp, "Atoms\n\n");
    int i = 0;
    for (unitedAtom ua : polymerChains)
        fprintf(fp, "%d\t%d\t%d\t%lf\t%lf\t%lf\t%d %d %d\n", ++i, ua.chainID, ua.type, ua.pos[0], ua.pos[1], ua.pos[2], 
		(int)(ua.pos[0]/boxSize[0]), (int)(ua.pos[1]/boxSize[1]), (int)(ua.pos[2]/boxSize[2]));
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
    for (int i = 0; i < nDihedrals(); i++)
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

void SARW::report() const
{
    logPrintf("\nSUMMARY:\n");
        
    logPrintf("No. of united atoms: ");
    logPrintf("desired = %d, ", targetCount());
    logPrintf("actual = %d\n", actualCount());
    
    logPrintf("Number density: ");
    logPrintf("desired = %.3e, ", targetNumberDensity());
    logPrintf("actual = %.3e\n", actualNumberDensity());

    logPrintf("Mass density: ");
    logPrintf("desired = %.2f, ", targetMassDensity);
    logPrintf("actual = %.2f\n", actualMassDensity());

    logPrintf("\nExporting distribution of chain lengths: chains.dat\n");
    FILE* fp = fopen("chains.dat", "w");
    for (int i : chainLengths) fprintf(fp, "%d\n", i);
    fclose(fp);
}

void SARW::setLogFlags(string logProgressFlag, string logStepsFlag)
{
    logProgress = logProgressFlag=="yes";
    logSteps = logStepsFlag=="yes";
    if (logSteps) logProgress=true;
}
