# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from scipy import io
import numpy as np
import pylab as pl
from scipy.signal import savgol_filter as sg


def smooth(x, fillen=27, polyord=3, axis=-1):
    y = sg(np.flip(sg(x, fillen, polyord, axis=axis), axis),
           fillen, polyord, axis=axis)
    return np.flip(y, axis)
    

dat = io.loadmat('TFRPowerSummary.mat')
f = dat['f'].flatten()
c7 = dat['c7']
c7 = c7.T - c7.mean(axis=1).T
c14 = dat['c14']
c14 = c14.T - c14.mean(axis=1).T

c20 = dat['c20']
c20 = c20.T - c20.mean(axis=1).T

delta = (smooth(c20, axis=0) - smooth(c14, axis=0)).T
dt = delta[:26, np.logical_and(f < 42, f > 13)].sum(axis=1)
dmt = dt.mean(axis=0)
da = delta[26:, np.logical_and(f < 42, f > 13)].sum(axis=1)
dma = da.mean(axis=0)
det = dt.std(axis=0) / (26 ** 0.5)
dea = da.std(axis=0) / (21 ** 0.5)

dt2 = delta[:26, np.logical_and(f > 42, f > 62)].sum(axis=1)
dmt2 = dt2.mean(axis=0)
da2 = delta[26:, np.logical_and(f > 42, f > 62)].sum(axis=1)
dma2 = da2.mean(axis=0)
det2 = dt2.std(axis=0) / (26 ** 0.5)
dea2 = da2.std(axis=0) / (21 ** 0.5)

# Plot spectrum
ntime = 50
delta = delta * ntime  # number of time bins pooled
dtf = np.mean(delta[:26, ], axis=0) 
daf = np.mean(delta[26:, ], axis=0) 

nsubj = delta.shape[0]
de = delta.std(axis=0) / (nsubj ** 0.5 * 2.)
cont = dtf - daf
pl.plot(f, cont, 'k')
pl.fill_between(f, cont + de, cont - de, color='k', alpha=0.2)
pl.xlabel('Frequency (Hz)', fontsize=14)
pl.ylabel('$\Delta P_{TD}$ - $\Delta P_{ASD}$ (dB)',
          fontsize=14)
pl.xlim((12, 65))
pl.ylim((-3., 3.))