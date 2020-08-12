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

if __name__ == "__main__":
	NA = 6.022e23
	L = np.array([42.,42.,336.])

	filename = 'polymer.dat'

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
