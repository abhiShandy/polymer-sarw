import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# v = np.loadtxt('chains.dat')

# len_min = min(v)
# len_max = max(v)

# plt.hist(v, facecolor='blue', alpha=1)

# plt.xlabel('No. of Chains')
# plt.ylabel('Chain Length')
# plt.title('Chain Length Distribution')
# # plt.axis([len_min, len_max, 1, len(v)])
# plt.grid(True)

# plt.savefig('dist.png')
# plt.savefig('dist.pdf')
# plt.show()
g = np.loadtxt('thermo.dat')

plt.plot(g[5:,3], g[5:,1])
plt.title('Glass Transition Temperature - Polyethylene')
plt.xlabel("Temperature (Kelvins)")
plt.ylabel("Volume (cubic Angstrom)")
plt.grid(True)


plt.savefig('tg.png')
plt.show()