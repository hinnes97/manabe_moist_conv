import numpy as np
import matplotlib.pyplot as plt

dat = np.loadtxt('moist_dat.txt')
P = dat[:,0]
Told = dat[:,1]
qold = dat[:,2]
true_ad = dat[:,3]
T_adj = dat[:,4]
true_ad_post = dat[:,5]
q_adj = dat[:,6]

dat2 = np.loadtxt('dry_dat.txt')
p_dry = dat2[:,0]
T_adj_dry = dat2[:,4]
true_ad_post2 = dat2[:,5]

plt.semilogy(Told,P, label='Initial profile')
plt.semilogy(T_adj, P, label= 'Adjusted profile')
plt.semilogy(true_ad, np.flip(P), label='Moist adiabat before adjustment')
plt.semilogy(true_ad_post, np.flip(P), label='Moist adiabat after adjustment')
plt.semilogy(T_adj_dry, p_dry, label='Adjusted using dry enthalpy conservation')
#plt.semilogy(true_ad_post2, np.flip(P), label='Moist adiabat after adjustment, dry')

plt.gca().invert_yaxis()
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Pa)')
plt.tight_layout()

plt.savefig('test_t.pdf')

plt.figure()

plt.loglog(qold, P, label='q before')
plt.loglog(q_adj, P, label='q after')
plt.xlabel('q (kg/kg)')
plt.ylabel('Pressure (Pa)')
plt.gca().invert_yaxis()
plt.legend()
plt.tight_layout()
plt.savefig('test_q.pdf')
