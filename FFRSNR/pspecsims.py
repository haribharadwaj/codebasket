# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 12:43:34 2013

@author: hari
"""
import numpy as np
from scipy.special import iv
# from scipy import floor,ceil
import matplotlib
matplotlib.use('PS')   # generate postscript output by default

import pylab as pl

def db2mag(x):
    """ Converts from dB to magnitute ratio
    
    Parameters
    ----------
    
    x - Input in dB
    
    Returns
    -------
    
    m - magnitude ratio
    """
    
    m = 10.0**(x/20.0)
    return m
    
def db(x):
    """ Converts to decibels
    
    Parameters
    ----------
    
    x - Input in linear units
    
    Returns
    -------
    
    y - Equivalend in decibel units
    """
    
    y = 20*np.log10(x)
    return y
    
def pspec(z,Npairs):
    """ Pairwise power estimate (single spectral bin)
    
    Parameters
    ----------
    
    z - Input array of shape (Ntrials,channels) or (Ntrials,)
    
    Npairs - Number of pairs of trials to estimate from
    
    Returns
    -------
    
    pS - Pairwise power estimate
    """
    ntrials = z.shape[0]
    trial_pairs = np.random.randint(0,ntrials,(Npairs,2))
    trial_pairs = trial_pairs[np.not_equal(trial_pairs[:,0],trial_pairs[:,1])]
    first = trial_pairs[:,0]
    second = trial_pairs[:,1]
    
    #first = np.arange(0,np.int(floor(ntrials/2.0)))
    #second = np.arange(np.int(ceil(ntrials/2.0)),ntrials)
    
    
    xw_1 = z[first]
    xw_2 = z[second]
    pS = np.real((xw_1*xw_2.conj()).mean(axis = 0))
    
    return pS
                     

def setcafontsize(fontsize = 20):
    """Set fontsize  for current axis
    
    Parameters
    ----------
    fontsize - Default is 20
    """
    ax = pl.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(20)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(20)
        
    
Ntrials = 100
Nsims = 5000
# Draw signal phase from von mises density
kappa = 1
phi_sig = np.random.vonmises(0, kappa, size = (Ntrials,Nsims))
PLV_true = iv(1,kappa)/iv(0,kappa)

phi_n = np.random.rand(Ntrials,Nsims)*2*np.pi


SNRlist = np.arange(-6.0,6.2,0.2) #dB

pS_mu = np.zeros(SNRlist.shape)
pS_std = np.zeros(SNRlist.shape)
S_mu = np.zeros(SNRlist.shape)
S_std = np.zeros(SNRlist.shape)
S = 10.0

for nSNR, SNR in enumerate(SNRlist):
    print 'Running simulations for SNR = ', SNR
        
    A = db2mag(SNR)
    
    # Measurement
    z = S*np.exp(1j*phi_sig) + (S/A)*np.exp(1j*phi_n)
    
    Npairs = 4000
    pS = pspec(z,Npairs)
    pS_mu[nSNR] = pS.mean()
    pS_std[nSNR] = pS.std()
    S_mu[nSNR] = (np.abs(z.mean(axis=0))**2).mean()
    S_std[nSNR] = (np.abs(z.mean(axis=0))**2).std()
    
trueVal = (S**2)*(PLV_true**2)

pl.subplot(2,1,1)
pl.plot(SNRlist,pS_mu*100/trueVal,'k-',linewidth = 2)
pl.hold(True)
pl.plot(SNRlist,S_mu*100/trueVal,'g-',linewidth = 2)
pl.plot(SNRlist,np.ones(SNRlist.shape)*100,'r-',linewidth = 2)
pl.xlabel('SNR (dB)',fontsize = 20)
pl.ylabel('Power Estimate',fontsize = 20)
pl.xlim((np.min(SNRlist),np.max(SNRlist)))
pl.legend(('Mean Pairwise Estimate','Mean Spectral Estimate','Ground Truth'))
setcafontsize()

pl.subplot(2,1,2)
pl.plot(SNRlist,pS_std*100/trueVal,'k',linewidth = 2)
pl.hold(True)
pl.plot(SNRlist,S_std*100/trueVal,'g',linewidth = 2)
pl.xlabel('SNR (dB)',fontsize = 20)
pl.ylabel('Estimator standard deviation',fontsize = 20)
pl.xlim((np.min(SNRlist),np.max(SNRlist)))
pl.legend(('Pairwise Estimate','Spectral Estimate'))
setcafontsize()
pl.show()




