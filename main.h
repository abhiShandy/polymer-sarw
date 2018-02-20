const int nChains = 200;
const int minChainLength = 100;
const double mass_density = 0.92;	// units =  g/cc

const double N_avogadro = 6.022 * pow(10,23);
const int nUnitedAtoms = nChains * minChainLength;
const double number_density = N_avogadro*mass_density/14.0;	// units = per cc
const int boxSize = nearbyint(pow(10,8)*pow(nUnitedAtoms/number_density, 0.3333));
//#define POLYSTYRENE
//#define POLYPROPYLENE
#define POLYETHYLENE
#define DEBUG false
#define LOG true

const double minDist = 2.0;
const int maxTrials = 100;

// #ifndef
// #define TRIALS
// iTrial = 0;
// 	while(iTrial++ < maxTrials)
// 	{
// #endif
