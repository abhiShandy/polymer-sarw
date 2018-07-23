import numpy as np
import matplotlib.pyplot as plt
from itertools import islice

def nearestDist(x):
    dx = x[:,None,:] - x[None,:,:]
    dx -= np.floor(0.5 + dx) #wrap to nearest image convention
    dxSq = np.sum((L*dx)**2, axis=-1)
    iUpper = np.triu_indices(N, 1)
    return np.sqrt(np.min(dxSq[iUpper]))

#Unnormalized rdf between two pairs of atom positions (each Nx3)
def rdfPair(atoms1, atoms2, rMin, dr):
    dist = (atoms1[:,None,:] - atoms2[None,:,:]) * (1./L)
    dist -= np.floor(0.5 + dist)   #Wrap to nearest image convention
    dist *= L
    distMag = np.sqrt(np.sum(dist**2, axis=-1)).flatten()
    rMax = np.sqrt(np.sum((0.5*L)**2))
    rBinEdges = np.arange(rMin, rMax, dr)
    return np.histogram(distMag, rBinEdges)

#Uniform rdf:
def rdfUniform(rMin, dr):
	S = np.ceil(L/dr).astype(int)
	dx = L/S
	grids = [ np.arange(S[iDim])*dx[iDim] for iDim in range(3) ]
	meshes = np.meshgrid(*grids, indexing='ij')
	mesh = np.array([m.flatten() for m in meshes]).T
	#Subsample reference positions:
	nSample = 2.
	grids = [ np.arange(nSample)*dx[iDim]/nSample for iDim in range(3) ]
	meshes = np.meshgrid(*grids, indexing='ij')
	offsets = np.array([m.flatten() for m in meshes]).T
	#Accumulate RDFs for each offset - mesh combinations
	rdf,binEdges = rdfPair(mesh, offsets, rMin,dr)
	return rdf*(1./(mesh.shape[0]*offsets.shape[0])),binEdges

#Calculate intra and inter rdf of atom in chains specified by [chainID, x,y,z]
def rdfChains(atoms, rMin, dr):
	chainIDs = np.unique(atoms[:,0])
	atomChains = [ atoms[np.where(atoms[:,0]==chainID)[0], 1:] for chainID in chainIDs ]
	#Create normalizing rdf:
	rdfDen, binEdges = rdfUniform(rMin, dr)
	rdfDen *= (atoms.shape[0]**2) #TODO verify
	#Accumulate intra and inter combinations:
	rdfIntra = np.zeros(rdfDen.shape)
	rdfInter = np.zeros(rdfDen.shape)
	for iChain1,atoms1 in enumerate(atomChains):
		for iChain2,atoms2 in enumerate(atomChains):
			if iChain1==iChain2: #intra contribution
				rdfIntra += rdfPair(atoms1,atoms2, rMin,dr)[0]
			elif iChain1>iChain2: #inter contribution
				rdfInter += rdfPair(atoms1,atoms2, rMin,dr)[0]*2. #to account for iChain1<iChain2
	#Normalize rdfs:
	rdfIntra /= rdfDen
	rdfInter /= rdfDen
	return rdfIntra, rdfInter, rdfDen, binEdges

'''
def rdf(atoms,rdftype):
    
    atoms is an array saved in the form of [chainID, x, y, z]
    rdftype selects whether desired rdf is interchain or intrachain
    
    cum_rdf = [] # Store the cumulative average rdf based on inter/intra
    for atom in atoms:   # Iterate through all atoms
        print atom
        refAtom = atom              # Set this current atom as the reference atom
        dist = atoms-refAtom        # Find the distance b/w refAtom and all atoms
        dist[...,1:] = dist[...,1:]/L
        dist[...,1:] -= np.floor(0.5 + dist[...,1:])
        dist[...,1:] *= L
        rdf = []                    # Stores the rdf for the current refAtom
        if rdftype == 'intra':
            dist = dist[np.where(dist[...,0]==0)]  # get the distances from refAtom if chainID value after dist calc = 0 for intra rdf
        elif rdftype == 'inter':
            dist = dist[np.where(dist[...,0]!=0)]  # get the distances from refAtom if chainID value after dist calc != for inter rdf
        norm = np.linalg.norm(dist,axis=1)         # Get scalar distance
        rMax = np.amax(norm)                       # pick an arbitrary maximum radius to perform rdf
        rMax += 50                   # add a buffer to make sure that all the atoms are accounted for
        radii = np.arange(0,rMax)    # create an array of radii to iterate over for the rdf
        dr = 1.                      # use this assigned range to determine the rdf
        for r in radii:
            #normalize = (0.92 * NA / (1e21)) * 4 * np.pi * r**2 * dr
            rdf.append((len(norm[abs(norm -r - dr/2.) < dr/2.])))#/normalize)  # store in rdf the atoms within desired range of rdf
        rdf=np.array(rdf)            # convert to np array
        cum_rdf.append(rdf)          # store current rdf in a list
    cum_rdf = np.array(cum_rdf)      # convert to np array
    # Make all of the arrays within cum_rdf the same length
    # Different lengths due to different chain lengths in the system
    b = np.zeros([len(cum_rdf),len(max(cum_rdf,key=lambda x: len(x)))])
    for i,j in enumerate(cum_rdf):
        b[i][0:len(j)] = j
    mean_rdf = np.mean(b,axis=0)     # take the average of all the intra/inter rdfs for the system
    radii = np.arange(0,len(mean_rdf))  # make radii the same length for plotting

    return mean_rdf,radii

def uniform(radii):
    distr = np.random.uniform(0.,1.,(12391,3))
    refAtom = distr[0]
    dist = distr-refAtom
    dist -= np.floor(0.5 + distr)
    dist *= L
    norm = np.linalg.norm(dist,axis=1)
    rMax = np.amax(norm)
    rMax += 50
    radii = np.arange(0,rMax)
    dr = 1.
    rdf = []
    for r in radii:
        rdf.append((len(norm[abs(norm -r - dr/2.) < dr/2.])))
    rdf = np.array(rdf)
    return rdf,radii
'''

if __name__ == "__main__":
    NA = 6.022e23
    L = np.array([42.,42.,336.])

    filename = 'polymer.data'

    with open(filename) as lines:
        # Stores Chain ID, Position x, Position y, Position z
        atoms=np.genfromtxt(islice(lines,24,12391),usecols=(1,3,4,5))

    # STORING DATA
    # Get chain lengths stored in the chains.dat file
    filename = 'chains.dat'
    chainLen=np.genfromtxt(filename)
    chainLen[0]-=1
    atoms = atoms[atoms[:,0].argsort()] # Sort atoms by chain ID
    
    rMin = 2. #ignore self and nearest neighbor
    dr = 0.5
    rdfIntra,rdfInter,wr,binEdges = rdfChains(atoms, rMin,dr)
    binCenters = 0.5*(binEdges[:-1]+binEdges[1:])
    plt.plot(binCenters, rdfIntra, 'r', label='Intra')
    plt.plot(binCenters, rdfInter, 'b', label='Inter')
    plt.legend()
    plt.xlabel('Radius (nm)')
    plt.ylabel('Fraction of total atoms in system')
    plt.xlim([0,None])
    plt.ylim([0,None])
    plt.show()

    ## PLAN OF ATTACK
    # Vectorize each of the chains
    # Find the self rdf of the chain
    # Find the rdf of all other atoms
    # N x N x Nbins (in this case radii) distances
    # Normalize by uniform distribution
    # Nideal = N/V * Vshell(r)
    # Nideal = rho * (4*pi*r**2) dr
