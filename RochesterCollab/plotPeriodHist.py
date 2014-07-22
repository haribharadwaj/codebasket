# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 09:46:47 2014

@author: hari
"""
import numpy as np
from scipy import io
import pylab as pl


pl.rcParams.update({'font.size': 14})
froot = '/home/hari/Documents/RochesterCollab/'
bf = '3500'

SNRs = map(str, [1, 2, 3])
actualSNR = [20, 35, 50]

for k, SNR in enumerate(SNRs):
    fname = 'Hari_' + bf + 'Hzcell_SNR' + SNR + 'data.mat'
    dat = io.loadmat(froot + fname)
    x = np.squeeze(dat['x']/1000.0)
    y = np.squeeze(dat['y'])
    pl.figure(num=1)
    pl.subplot(3, 1, k+1)
    pl.plot(x, y, '.k', markersize=8)
    pl.xlim((0, 500))
    pl.ylabel('Trial #')
    pl.title('SNR = ' + str(actualSNR[k]) + ', BF = ' + bf)

    x = x[np.logical_and(x > 25, x < 500)]
    nfolds = np.round((500-25)/10)
    x_fold = np.remainder(x, 10.0)
    nbins = 25
    count, edges = np.histogram(x_fold, bins=nbins)
    binsize = np.mean(np.diff(edges))/1000.0
    cyc = (edges[1:] + edges[:-1])*0.5/10.0
    pl.figure(num=2)
    pl.subplot(3, 1, k+1)
    ntrials = 100
    pl.plot(cyc, count / (binsize * ntrials * nfolds), 'ko-', linewidth=2)
    pl.ylim((2, 150))
    pl.ylabel('Rate')
    pl.title('SNR = ' + str(actualSNR[k]) + ', BF = ' + bf)
    vs = np.abs(np.sum(count*np.exp(1j*cyc*2*np.pi)))/count.sum()
    vecstr = 'Vector strength = ' + str(np.round(vs, decimals=3))
    pl.text(0.01, 70, vecstr)

pl.figure(num=1)
pl.xlabel('PST (ms)')
pl.show()

pl.figure(num=2)
pl.xlabel('Period')
pl.show()
