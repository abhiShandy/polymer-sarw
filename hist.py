import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

v = np.loadtxt('chains.dat')

len_min = min(v)
len_max = max(v)

plt.hist(v, facecolor='blue', alpha=1)
# plt.hist(v, 50, normed=1, facecolor='blue', alpha=0.75)

plt.xlabel('No. of Chains')
plt.ylabel('Chain Length')
plt.title('Chain Length Distribution')
# plt.axis([len_min, len_max, 1, len(v)])
plt.grid(True)

plt.savefig('dist.png')
plt.savefig('dist.pdf')
plt.show()