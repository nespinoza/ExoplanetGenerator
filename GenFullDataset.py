# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import numpy as np
import batman

def init_batman(t,law):
    """
    This function initializes the batman code.
    """
    params = batman.TransitParams()
    params.t0 = 0.
    params.per = 1.
    params.rp = 0.1
    params.a = 15.
    params.inc = 87.
    params.ecc = 0.
    params.w = 90.
    params.u = [0.1,0.3]
    params.limb_dark = law
    m = batman.TransitModel(params,t)
    return params,m

def get_transit_model(t,t0,P,p,a,inc,q1,q2,ld_law):
    params,m = init_batman(t,law=ld_law)
    coeff1,coeff2 = reverse_ld_coeffs(ld_law, q1, q2)
    params.t0 = t0
    params.per = P
    params.rp = p
    params.a = a
    params.inc = inc
    params.u = [coeff1,coeff2]
    return m.light_curve(params)

def convert_ld_coeffs(ld_law, coeff1, coeff2):
    if ld_law == 'quadratic':
        q1 = (coeff1 + coeff2)**2
        q2 = coeff1/(2.*(coeff1+coeff2))
    elif ld_law=='squareroot':
        q1 = (coeff1 + coeff2)**2
        q2 = coeff2/(2.*(coeff1+coeff2))
    elif ld_law=='logarithmic':
        q1 = (1-coeff2)**2
        q2 = (1.-coeff1)/(1.-coeff2)
    return q1,q2

def reverse_ld_coeffs(ld_law, q1, q2):
    if ld_law == 'quadratic':
        coeff1 = 2.*np.sqrt(q1)*q2
        coeff2 = np.sqrt(q1)*(1.-2.*q2)
    elif ld_law=='squareroot':
        coeff1 = np.sqrt(q1)*(1.-2.*q2)
        coeff2 = 2.*np.sqrt(q1)*q2
    elif ld_law=='logarithmic':
        coeff1 = 1.-np.sqrt(q1)*q2
        coeff2 = 1.-np.sqrt(q1)
    return coeff1,coeff2

import ajplanet as rv_model

t0 = 2457903.33611
P = 7.5
# Cycles of the transit datasets:
tr_ncycle = np.array([5,30])
tr_ncycles = np.array([4,0])
tr_phases = np.array([0.5,0.02])
tr_instrument_names = ['K2','LCOGT']
# Cycles of the RV datasets:
rv_ncycle = np.array([20,25])
rv_ncycles = np.array([1,1.5])
rv_instrument_names = ['HARPS','CORALIE']
rv_mu = np.array([1.5,10.0])
rv_texp = [2.0,1.0]
p = np.array([0.05,0.05])
a = np.array([13.0,13.0])
inc = 89.5
ecc = 0.3
omega = 85.0
# Semi-amplitude (m/s):
K = 80.0
# LD coeffs of lightcurves:
q1 = np.array([0.2,0.6])
q2 = np.array([0.4,0.2])
# Lightcurve noise (ppm):
sigma_w = np.array([200.,1000.])
# RV noise (m/s):
sigma_w_rv = np.array([5.0,30.0])
sampled = [True,False]
# Kepler sampling time and exposures per sampling, 
# see Murphy (2012), Gilliland et al. (2010b).
texp = [0.01881944,0.00347222]
Nsamp = [271,1] 
ld_law = 'quadratic'
fname_lc = 'example_lc.dat'
fname_rv = 'example_rvs.dat'

# Generate dataset. First, lightcurve:
fout_lc = open(fname_lc,'w')
for i in range(len(tr_ncycle)):
    if sampled[i]:
        tstart = t0 + (tr_ncycle[i]-tr_phases[i])*P 
        tend = t0 + (tr_ncycle[i]+tr_ncycles[i]+tr_phases[i])*P
        tcurrent = tstart
        t = np.array([])
        f = np.array([])
        while True:
            tsampling = np.array([])
            fsampling = np.array([])
            for k in range(Nsamp[i]):
                tcurrent = tcurrent + (texp[i]/Nsamp[i])
                tsampling = np.append(tsampling,tcurrent)
            fsampling = get_transit_model(tsampling,\
                        t0,P,p[i],a[i],inc,q1[i],q2[i],ld_law)
            t = np.append(t,np.mean(tsampling))
            f = np.append(f,np.mean(fsampling))
            if tcurrent >= tend:
                break
        f = f + np.random.normal(0,sigma_w[i]*1e-6,len(t))
    else:
        tstart = t0 + (tr_ncycle[i]-tr_phases[i])*P
        tend = t0 + (tr_ncycle[i]+tr_ncycles[i]+tr_phases[i])*P
        t = np.arange(tstart,tend,texp[i])
        fsampling = get_transit_model(t,\
                    t0,P,p[i],a[i],inc,q1[i],q2[i],ld_law)
        f = fsampling + np.random.normal(0,sigma_w[i]*1e-6,len(t))
    for ii in range(len(t)):
            fout_lc.write('{0:.10f} {1:.10f} {2:.10f} '.format(t[ii],f[ii],sigma_w[i]*1e-6)+tr_instrument_names[i]+'\n')

fout_lc.close()
fout_rv = open(fname_rv,'w')
# Now RV data:
for i in range(len(rv_ncycle)):
        tstart = t0 + (rv_ncycle[i])*P
        tend = t0 + (rv_ncycle[i]+rv_ncycles[i])*P
        t = np.arange(tstart,tend,rv_texp[i]) 
        t = t + np.random.uniform(-0.25,0.25,len(t))
        rvs = rv_model.pl_rv_array(t,rv_mu[i],K,\
                                omega*np.pi/180.,ecc,\
                                t0,P)
        rvs = rvs + np.random.normal(0,sigma_w_rv[i],len(t))
        for ii in range(len(t)):
            fout_rv.write('{0:.10f} {1:.10f} {2:.10f} '.format(t[ii],rvs[ii],sigma_w_rv[i])+rv_instrument_names[i]+'\n')
fout_rv.close()
