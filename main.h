const int nChains = 20;
const int minChainLength = 1000;
const double mass_density = 0.92;	// units =  g/cc

const double N_avogadro = 6.022 * pow(10,23);
const int nUnitedAtoms = nChains * minChainLength;
const double number_density = N_avogadro*mass_density/14.0;	// units = per cc
int boxSize = nearbyint(pow(10,8)*pow(nUnitedAtoms/number_density, 0.3333));


// const int boxSize = 80;
//#define POLYSTYRENE
//#define POLYPROPYLENE
#define POLYETHYLENE
#define DEBUG false
#define LOG true
// #define CENTER_SEED

const int seed_radius = 10;

const double minDist = 2.0;
const int maxTrials = 100;

// #ifndef
// #define TRIALS
// iTrial = 0;
// 	while(iTrial++ < maxTrials)
// 	{
// #endif
