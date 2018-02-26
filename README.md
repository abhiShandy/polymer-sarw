## Description of Model Parameters

All model parameters are in the main.h file

- nChains = number of chains. These many number of seeds are randomly placed in the allowable region of simulation boxx

- minChainLength = target length of each chain. Propagation of chains do not stop when they reach this length.

- mass_density = desired mass density of polymer (g/cc)

- minDist = minimum distance defined between non-bonded monomers

- maxTrials = upper limit on the MC steps