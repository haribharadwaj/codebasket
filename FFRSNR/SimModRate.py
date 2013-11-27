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
datdir = '/home/hari/Documents/PythonCodes/research/FFRSNRstuff/'
x = io.loadmat(datdir + 'x.mat')['x'][30,:,:]

fs = 4096.0
ntime = x.shape[1]
fsig = 211.0 #Picking some frequency nowhere neaR 100 Hz
Strue = 0.1e-6 #Typical FFR size
t = np.arange(0.0,ntime)/fs
sig = Strue*np.sin(2*np.pi*fsig*t)

noisefac_list = np.arange(0.2,10.0,0.2)
SIlist = np.zeros(noisefac_list.shape)
Srawlist = np.zeros(noisefac_list.shape)
pspeclist = np.zeros(noisefac_list.shape)
for k, noisefac in enumerate(noisefac_list):
    print 'Testing Noise Level #',k
    y = sig + x*noisefac
    # Comuting spectrum and phase-locking
    params = dict(Fs = fs, tapers = [1, 1], pad = 1,
              fpass = [5,400], Npairs = 5000, itc = 1)
    [Sraw,f] = spectral.mtspecraw(y, params)
    [pspec,f] = spectral.mtpspec(y, params)
    Srawlist[k] = max(Sraw[abs(f - fsig) < 1])
    pspeclist[k] = max(pspec[abs(f - fsig) < 1]) 
    
    
pl.plot(noisefac_list,Srawlist,linewidth = 2)
pl.figure()
pl.plot(noisefac_list,pspeclist,'k',linewidth = 2)
pl.show()