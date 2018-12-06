#ifndef SARW_PERIODIC_LOOKUP_H
#define SARW_PERIODIC_LOOKUP_H

inline vector3<> mul(const vector3<>& a, const vector3<>& b)
{	return vector3<>(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}

//! O(1) look-up table for finding periodic image within distance
class PeriodicLookup
{	vector3<> boxSize; double thresh;
	vector3<> boxSizeInv;
	vector3<int> S; //lookup mesh sample count
	std::vector<vector3<>> points;
	std::vector< std::vector<size_t> > indices; //list of indices into points, for each lookup mesh cell
	
	inline size_t meshIndex(vector3<int> iv) const
	{	for(int k=0; k<3; k++) //wrap to [0,S)
		{	if(iv[k] < 0) iv[k] += S[k];
			if(iv[k] >= S[k]) iv[k] -= S[k];
		}
		return iv[0]+S[0]*size_t(iv[1]+S[1]*iv[2]);
	}

public:
	//!Initialize given box size and a threshold for distance detection
	PeriodicLookup(vector3<> boxSize, double thresh) : boxSize(boxSize), thresh(thresh)
	{	//Set S such that bin size just exceeds thresh in each dimensions:
		for(int k=0; k<3; k++)
		{	S[k] = int(floor(boxSize[k]/thresh));
			boxSizeInv[k] = 1./boxSize[k];
		}
		//Initialize indices:
		indices.resize(S[0]*S[1]*S[2]);
	}
	
	//! Add point updating indices
	void addPoint(const vector3<> pos)
	{	vector3<> v = mul(pos, boxSizeInv); //to lattice coordinates
		vector3<int> iv;
		for(int k=0; k<3; k++)
		{	v[k] -= floor(v[k]); //in [0,1)
			iv[k] = int(floor(v[k]*S[k] + 0.5)); //in [0, S]
		}
		indices[meshIndex(iv)].push_back(points.size());
		points.push_back(v);
	}
	
	//! Return whether there are any atoms within thresh of the specified location,
	//! optionally ignoring atom ignoreIndex if >= 0
	bool find(vector3<> pos, int ignoreIndex=-1) const
	{	vector3<int> ivMin, ivMax;
		vector3<> v = mul(pos, boxSizeInv); //to lattice coordinates
		for(int k=0; k<3; k++)
		{	v[k] -= floor(v[k]); //in [0,1)
			double thresh_k = thresh * boxSizeInv[k]; //max spacing along this direction in lattice coordinates
			ivMin[k] = int(floor((v[k]-thresh_k)*S[k] + 0.5)); //in [-1, S]
			ivMax[k] = int(floor((v[k]+thresh_k)*S[k] + 0.5)); //in [0, S]
		}
		vector3<int> iv;
		for(iv[0]=ivMin[0]; iv[0]<=ivMax[0]; iv[0]++)
		for(iv[1]=ivMin[1]; iv[1]<=ivMax[1]; iv[1]++)
		for(iv[2]=ivMin[2]; iv[2]<=ivMax[2]; iv[2]++)
			for(size_t index: indices[meshIndex(iv)])
				if(ignoreIndex<0 || int(index)!=ignoreIndex)
				{	//Compute distance:
					vector3<> dv = v - points[index];
					for(int k=0; k<3; k++)
						dv[k] -= floor(0.5 + dv[k]); //wrap to [-0.5,0.5); minimum image convention
					double dist = mul(boxSize, dv).length();
					if(dist <= thresh)
						return true; // found neighbour within threshold
				}
		return false; //no neighbours found
	}
};

#endif //SARW_PERIODIC_LOOKUP_H