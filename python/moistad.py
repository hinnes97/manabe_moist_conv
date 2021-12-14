import numpy as np
import matplotlib.pyplot as plt
import phys
from scipy import integrate

def satvp(T,wet):
    """Calculates saturation vapour pressure"""
    T0 = wet.TriplePointT
    p0 = wet.TriplePointP

    if T>wet.TriplePointT:
        L = wet.L_vaporization
    else:
        L = wet.L_sublimation
        
    return p0*np.exp( -L/wet.R*(1/T - 1/T0) )
    
def dlntdlnp(lnp,lnT,dry,wet):
    """Calculates moist adiabatic gradient (see Ding 2016)
       This value is correct, even in the nondilute regime
    """
    T = np.exp(lnT)
    p = np.exp(lnp)
    
    dry.update()
    wet.update()

    Rd,Rcpd,cpd = dry.R,dry.Rcp,dry.cp
    Rw,cpw = wet.R,wet.cp
    eps = wet.MolecularWeight/dry.MolecularWeight

    if T>wet.TriplePointT:
        L = wet.L_vaporization
    else:
        L = wet.L_sublimation

    psat = satvp(T,wet)
    qsat = eps*psat/p
    pa = p - psat
    
    rsat = qsat/(1-qsat)
    
    num = 1 + (L/Rd/T)*rsat
    den = 1 + ( (cpw/cpd) + ( (L/Rw/T) - 1 )*(L/cpd/T) )*rsat
    
    F = Rcpd*num/den

    dlnpdlnt = (psat/p)*L/Rw/T + pa/p/F
    
    return 1./dlnpdlnt

def moistadiabat(p,Ts,dry,wet):
    """Numerically integrates a moist adiabat from surface temp and pressure"""
    ps = p[0]
    pend = p[-1]
    lnpend = np.log(pend)

    lnps,lnTs = np.log(ps),np.log(Ts)

    sol = integrate.solve_ivp(dlntdlnp,(lnps,lnpend),[lnTs],args=(dry,wet),t_eval=np.log(p))

    pad = np.exp(sol.t)
    Tad = np.exp(sol.y)[0]
    return Tad,pad

