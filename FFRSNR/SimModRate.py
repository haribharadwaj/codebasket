# -*- coding: utf-8 -*-
"""
Simulating the adjusted modulated rate metric
Created on Sat Nov 23 14:09:28 2013

@author: Hari Bharadwaj
"""

from anlffr import spectral
import numpy as np
from scipy import io
import pylab as pl


# Load data recorded with a 100 Hz transposed tone
# datdir = '/home/hari/Documents/PythonCodes/research/FFRSNRstuff/'
datdir = '/home/hari/Documents/MNEforBiosemi/'
x = io.loadmat(datdir + 'x.mat')['x'][30,:,:]

fs = 4096.0
ntime = x.shape[1]
fsig = 211.0 #Picking some frequency nowhere neaR 100 Hz
Strue = 0.1e-7 #Typical FFR size
t = np.arange(0.0,ntime)/fs
sig = Strue*np.sin(2*np.pi*fsig*t)

noisefac_list = np.arange(0.2,5.0,0.2)
Slist = np.zeros(noisefac_list.shape)
pspeclist = np.zeros(noisefac_list.shape)
for k, noisefac in enumerate(noisefac_list):
    print 'Testing Noise Level #',k
    y = sig + x*noisefac
    # Comuting spectrum and phase-locking
    params = dict(Fs = fs, tapers = [2, 3], pad = 1,
              fpass = [5,400], Npairs = 1000, itc = 1)
    [S,N, f] = spectral.mtspec(y, params)
    [pspec,vpspec,f] = spectral.bootfunc(y, 100, 50, params, func = 'pspec')
    Slist[k] = max(S[abs(f - fsig) < 1])
    pspeclist[k] = max(pspec[abs(f - fsig) < 1]) 
    
    
pl.plot(noisefac_list,Slist,'b',linewidth = 2)
pl.hold(True)
pl.plot(noisefac_list,pspeclist**0.5,'r',linewidth = 2)
pl.ylim(0,1e-4)
pl.show()