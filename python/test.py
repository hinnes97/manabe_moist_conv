import numpy as np
import moistad as mad
import adjust as adj
import matplotlib.pyplot as plt
import phys

P = np.logspace(0,5,100)
dP = np.r_[np.diff(P)[0], np.diff(P)]
T = np.linspace(-600,300,100)
T[T<150]= 150
q = np.ones(100)*0.1

# Store old versions of T,q
Told = np.copy(T)
qold = np.copy(q)

true_ad, P_ad = mad.moistadiabat(np.flip(P), T[-1], phys.N2, phys.H2O)

T_adj, q_adj, mask = adj.moistAdj(T, P, dP, q, phys.N2, phys.H2O)

true_ad_post, P_ad_post = mad.moistadiabat(np.flip(P), T_adj[-1], phys.N2, phys.H2O)

np.savetxt('moist_dat.txt',
           [[P[i], Told[i], qold[i], true_ad[i], T_adj[i], true_ad_post[i], q_adj[i]] for i in range(len(T))])

#plt.savetxt('T_moist.txt', T_adj)
print('DONE', mask)
