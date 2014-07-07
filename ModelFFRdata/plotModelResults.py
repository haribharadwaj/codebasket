# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 12:35:30 2014

@author: hari
"""
import numpy as np
from scipy import io
import pylab as pl

subjlist = ['I13', 'I33', 'I41', 'I08', 'I02']

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/ModelData/'

pl.figure(num=1)
pl.figure(num=2)
for subj in subjlist:
    respath = froot + subj + '/RES/'
    # Magnitude plots
    condlist = np.arange(1, 16)
    f = (condlist-1)*30 + 100
    condstemlist = np.array(map(str, f))
    plv = np.zeros(condlist.shape)
    P = np.zeros(condlist.shape)
    for k, cond in enumerate(condstemlist):
        fname = respath + subj + '_' + cond + 'Hz_results.mat'
        f0 = f[k]
        dat = io.loadmat(fname)
        f_fft = dat['f'].squeeze()
        P[k] = 10 * (np.log10(dat['cpow'].squeeze()
                     [np.argmin(abs(f_fft-f0))]/1e-12))
        plv[k] = dat['cplv'].squeeze()[np.argmin(abs(f_fft-f0))]

    pl.figure(num=1)
    ax1 = pl.subplot(2, 1, 1)
    pl.hold(True)
    pl.plot(f, P, linewidth=2)
    pl.ylabel('Power (dB re: 1uV)', fontsize=16)
    pl.subplot(2, 1, 2, sharex=ax1)
    pl.hold(True)
    pl.plot(f, plv, linewidth=2)
    pl.ylabel('PLV (squared)', fontsize=16)
    pl.xlabel('Modulation frequency (Hz)', fontsize=16)

    # Phase plots
    fname = respath + subj + '_all_phase.mat'
    dat = io.loadmat(fname)
    ph = dat['Ph_f0'].squeeze()
    # phclean = np.median(unwrap(ph, axis=1), axis=0) / (2 * np.pi)
    phclean = np.unwrap(ph[30, :]) / (2 * np.pi)
    f = dat['f0_list'].squeeze()
    t_grp = -1e3 * np.diff(phclean) / (np.diff(f).mean())
    f_grp = f[1:] - (np.diff(f).mean())/2

    pl.figure(num=2)
    ax2 = pl.subplot(2, 1, 1)
    pl.hold(True)
    pl.plot(f, phclean, linewidth=2)
    pl.ylabel('Unwrapped phase (cycles)', fontsize=16)
    pl.subplot(2, 1, 2, sharex=ax2)
    pl.plot(f_grp, t_grp, linewidth=2)
    pl.hold(True)
    pl.xlabel('Modulation frequency (Hz)', fontsize=16)
    pl.ylabel('Group Delay (ms)', fontsize=16)
    pl.xlim(np.min(f), np.max(f))
    pl.legend(subjlist, loc=1)
pl.show()
