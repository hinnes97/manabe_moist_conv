import numpy as np
import matplotlib.pyplot as plt

py_dat = np.loadtxt('../python/moist_dat.txt')

p = py_dat[:,0]
T_adj = py_dat[:,4]

for_dat = np.loadtxt('output.txt')

p_for = for_dat[:,0]
Told  = for_dat[:,1]
T_adj_for = for_dat[:,2]

plt.semilogy( Told, p_for, label='Before')
plt.semilogy(T_adj_for,p_for, label='After')
plt.semilogy( T_adj, p, label='python')

plt.gca().invert_yaxis()
plt.xlabel('T (K)')
plt.ylabel('P (Pa)')
plt.tight_layout()

plt.savefig('fortran.pdf')
