# -*- coding: utf-8 -*-
"""
Simulating the adjusted modulated rate metric
Created on Sat Nov 23 14:09:28 2013

@author: Hari Bharadwaj
"""

from anlffr import spectral
import numpy as np
import pylab as pl

fs = 4096.0
T = 0.5
ntime = 4096.0*T
fsig = 211.0
Strue = 0.2
t = np.arange(0.0,ntime)/fs
sig = Strue*np.sin(2*np.pi*fsig*t)
Ntrials = 500
x = np.random.randn(Ntrials,ntime)
noisefac_list = np.arange(0.5,10.0,0.5)
Slist = np.zeros(noisefac_list.shape)
pspeclist = np.zeros(noisefac_list.shape)
for k, noisefac in enumerate(noisefac_list):
    
    print 'Testing Noise Level #',k
    y = sig + x*noisefac
    # Comuting spectrum and phase-locking
    TW = 1
    Ntaps = 1
    params = dict(Fs = fs, tapers = [TW, Ntaps], pad = 1,
              fpass = [5,400], Npairs = 5000, itc = 1)
    [S,N,f] = spectral.mtspec(y, params)
    [pspec,f] = spectral.mtpspec(y, params)
    
    bw = 2*TW/T
    Slist[k] = (S[abs(f - fsig) < bw]).mean()
    pspeclist[k] = (pspec[abs(f - fsig) < bw]).mean() 
    pl.figure(1)
    pl.plot(f, pspec)
    pl.hold(True)
    
    pl.figure(2)
    pl.plot(f, S)
    pl.hold(True)
        
pl.figure()    
pl.plot(noisefac_list,Slist,'b',linewidth = 2)
pl.hold(True)
pl.plot(noisefac_list,pspeclist**0.5,'r',linewidth = 2)
pl.ylim(1,25)
pl.show()