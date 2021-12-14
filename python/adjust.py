'''
HII December 2021 - neatened version of dilute moist adiabatic adjustment
'''

import numpy as np
import phys
import moistad
from scipy.optimize import root_scalar

def cond(T,p, dry, wet):
    '''Relax q to saturation vapour concentration'''
    eps = wet.MolecularWeight/dry.MolecularWeight
    
    Lvap_R = wet.L_vaporization_TriplePoint*wet.MolecularWeight/phys.Rstar
    Lsub_R = wet.L_sublimation*wet.MolecularWeight/phys.Rstar
    Tref = wet.TriplePointT
    pref = wet.TriplePointP

    qsat = lambda L_R : pref*np.exp(-L_R*(1/T - 1/Tref))/p * eps

    if T>Tref:
        qsats = qsat(Lvap_R)
    else:
        qsats = qsat(Lsub_R)

    return qsats
    #q = np.minimum(q0, qsats)

    # Cold trap ensures q above cold point doesn't exceed value at coldest point
    #index = signal.argrelmin(T,order=2)[0][0]
    #q[:index] = q[index]
    
    #return q
    

def moistAdj(T, p, dp, q, dry, wet):
    """Performs moist adjustment, conserving moist enthalpy
    The value of the gradient dlnT/dlnp is calculated in moistad.py
    We check temperatures pairwise going downwards and upwards through the atmosphere
    iteratively until convergence is reached
    We relate old temps (T1, T2) to new temps (T1', T2') by conserving enthalpy:
    T1*dp1 + T2*dp2 + L*q1*dp1 + L*q2*dp2 = T1'*dp1 + T2'*dp2 + L*q1'*dp1 + L*q2'*dp2
                                           = T1'(dp1 + (p2/p1)**(dlnT/dlnp)*dp2) + L*q1'*dp1 + L*q2'*dp2
    Therefore: 
    T1' = (T1*dp1 + T2*dp2 + L/cp*(q1-q1')*dp1 + L/cp*(q2-q2')*dp2)/(dp1 + (p2/p1)**(dlnT/dlnp)*dp2)
    T2' = T1'*(p2/p1)**(dlnT/dlnp)

    Note this has made the assumption that moist component of the enthalpy is 
    negligible (i.e. Lq terms are neglected). If we include these terms, a Newton
    iteration would have to be performed since new values of q, q', are dependent 
    on the temperatures T1' and T2' via Clausius Clapeyron. 
    
    For more info on non-dilute convection, see Ding 2016
    """
    #Leconte 2018 composition effect
    #vap = phys.H2O
    #vap.update()
    #omega = 1 - phys.H2.MolecularWeight/vap.MolecularWeight

    # Cold trap
    #qsats = np.array([cond(T,100*p) for (p,T) in zip(atm.p,atm.temp)])
    #mins = signal.argrelmin(qsats,order=20)[0]
    #if mins.size>0:
    #    indices = mins[q0>qsats[mins]]
    #    if indices.size>0:
    #        index=indices[-1]
    #        atm.mixing_ratios[0][:index+1] = qsats[index]
    #    else:
    #        index = 0
    #else:
    #    index = 0
    
    conv_mask = np.zeros_like(T)
    # Delta factor speeds up convergene
    delta = 0.0000001
    # play around with this number, trade-off between speed and convergence
    N_iter = 100

    for n in range(N_iter):
        print('n = ', n)
        # Downwards pass
        #print('MOIST ENTHALPY1:', moist_enthalpy(T, p, q, dp, dry, wet))
        for i in np.arange(len(T)-1):
            qsat1 = cond(T[i], p[i], dry, wet) 
            qsat2 = cond(T[i+1], p[i+1], dry, wet)
            
            if (q[i] > qsat1*(1-delta)) and (q[i+1] > qsat2*(1-delta)):
                q[i] = qsat1
                q[i+1] = qsat2
                T1,p1,dp1 = T[i],p[i],dp[i]
                T2,p2,dp2 = T[i+1],p[i+1],dp[i+1]
                pfact = (p1/p2)**moistad.dlntdlnp(np.log(p2), np.log(T2), dry, wet) 

                if T2 > wet.TriplePointT:
                    L = wet.L_vaporization
                else:
                    L = wet.L_sublimation

                if T1 < T2*pfact*(1+delta):
                    #Tbar = (dp1*T1+dp2*T2)/(dp1+dp2) 
                    #T2 = (dp1+dp2)*Tbar/(dp2+dp1*pfact)
                    #T1 = T2*pfact
                    T2 = root_scalar(iterator, args=(T1, T2, p1, p2, dp1, dp2, qsat1, qsat2, pfact, dry, wet), x0=250,x1=300).root
                    T1 = T2*pfact
                    T[i] = T1
                    T[i+1] = T2
                    q[i]   = cond(T1,p1, dry, wet)
                    q[i+1] = cond(T2,p2, dry, wet)
                    conv_mask[i] = 1
                    conv_mask[i+1] = 1
        #Upward pass
        for i in np.arange(len(T)-2,-1,-1):
            qsat1 = cond(T[i],p[i], dry, wet)
            qsat2 = cond(T[i+1],p[i+1], dry, wet)
            if (q[i]>qsat1*(1-delta)) and (q[i+1] > qsat2*(1-delta)):
                q[i] = qsat1
                q[i+1] = qsat2
                T1,p1,dp1 = T[i],p[i],dp[i]
                T2,p2,dp2 = T[i+1],p[i+1],dp[i+1]
                pfact = (p1/p2)**moistad.dlntdlnp(np.log(p2), np.log(T2),dry,wet)
                
                if T2 > wet.TriplePointT:
                    L = wet.L_vaporization
                else:
                    L = wet.L_sublimation
                    
                if T1 < T2*pfact*(1+delta):
                    #Tbar = (dp1*T1+dp2*T2)/(dp1+dp2)
                    #T2 = (dp1+dp2)*Tbar/(dp2+dp1*pfact)
                    #T1 = T2*pfact
                    
                    T2 = root_scalar(iterator, args=(T1, T2, p1, p2, dp1, dp2, qsat1, qsat2, pfact, dry, wet), x0=250,x1=300).root
                    T1 = T2*pfact

                    #T2 = (dp1+dp2)*Tbar/(dp2+dp1*pfact)
                    #T1 = T2*pfact
                    T[i] = T1
                    T[i+1] = T2
                    q[i]   = cond(T1,p1, dry, wet)
                    q[i+1] = cond(T2,p2, dry, wet)
                    conv_mask[i] = 1
                    conv_mask[i+1] = 1

        
    return T,q, conv_mask

def iterator(T2_new, T1, T2, p1, p2, dp1, dp2, q1, q2, pfact, dry, wet):
    cp = dry.cp
    
    if T2 > wet.TriplePointT:
        L = wet.L_vaporization
    else:
        L = wet.L_sublimation

    q1_new = cond(T2_new*pfact, p1, dry, wet)
    q2_new = cond(T2_new, p2, dry, wet)

    # Set L=0 to compare with convection scheme conserving dry enthalpy
    #L=0
    return 1 - (T1*dp1 + T2*dp2 + L/cp*(q1-q1_new)*dp1 + L/cp*(q2-q2_new)*dp2)/(dp2 + dp1*pfact) / T2_new

def moist_enthalpy(T, p, q, dp, dry, wet):
    L = np.zeros_like(T)
    cp = dry.cp
    L[T>wet.TriplePointT] = wet.L_vaporization
    L[T<wet.TriplePointT] = wet.L_sublimation

    # Set L=0 to compare with convection scheme conserving dry enthalpy
    #L=0
    
    return ((cp*T + L*q)*dp).sum()

