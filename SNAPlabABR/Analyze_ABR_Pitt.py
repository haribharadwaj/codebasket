#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 14:13:16 2024

@author: hari
"""

from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch
from scipy.signal import savgol_filter as sg
import pylab as pl
from scipy.io import savemat


# Setup bayesian-weighted averaging
def bayesave(x, trialdim=0, timedim=1, method='mean', smoothtrials=19):
    ntime = x.shape[timedim]
    wts = 1 / np.var(x, axis=timedim, keepdims=True)
    wts = sg(wts, smoothtrials, 3, axis=trialdim)  # Smooth the variance
    normFactor = wts.sum()
    wts = wts.repeat(ntime, axis=timedim)
    ave = (x * wts).sum(axis=trialdim) / normFactor
    return ave


# Adding Files and locations
# froot = '/Users/hari/Dropbox/Data_Testing/SpeechCoding/Data/ABR_EFR/'
froot = '/Users/HMB105/Library/CloudStorage/Dropbox/Data_Testing/SpeechCoding/Data/ABR_EFR/'


subj = 'HL011'
fpath = froot

ear = False  # True for left, False for right

if ear:
    cond = [1, 3]
    condname = 'Left'
    refchans = ['EXG1',]
    savestem = '_L_loud'
else:
    cond = [2, 4]
    condname = 'Right'
    refchans = ['EXG2',]
    savestem = '_R_loud'

print(f'Running Subject {subj}')

bdfs = fnmatch.filter(os.listdir(fpath), subj + '_ABR*.bdf')

if len(bdfs) >= 1:
    for k, bdf in enumerate(bdfs):
        edfname = fpath + bdf
        # Load data and read event channel
        raw, eves = bs.importbdf(edfname, refchans=refchans)

        # Pick channels to not include in epoch rejection
        raw.info['bads'] += ['EXG1', 'EXG2', 'CP2']
        # Filter the data
        raw.filter(l_freq=30., h_freq=3000)

        # Epoch the data
        tmin, tmax = -0.002, 0.080
        rejthresh = 200e-6  # Because of high-pass but using median
        epochs = mne.Epochs(raw, eves, cond, tmin=tmin, proj=False,
                            tmax=tmax, baseline=(-0.002, 0.0016),
                            reject=dict(eeg=rejthresh),
                            verbose='WARNING')
        xtemp = epochs.get_data()
        t = epochs.times * 1e3 - 1.6 # Adjust for delay and use ms
        # Reshaping so that channels is first
        if(xtemp.shape[0] > 0):
            xtemp = xtemp.transpose((1, 0, 2))
            if(k == 0):
                x = xtemp
            else:
                x = np.concatenate((x, xtemp), axis=1)
        else:
            continue
else:
    RuntimeError('No BDF files found!!')

    

# Purdue channels
# goods = [1, 5, 6, 27, 28, 30, 31]

# Pitt channels
goods = [28, 3, 30, 26, 4, 25, 7, 31, 22, 9, 8, 21, 11, 12, 18]


y = bayesave(x[goods, :, :].mean(axis=0)) * 1e6  # Microvolts


pl.plot(t, y, linewidth=2)
pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel(r'Response ($\mu$V)', fontsize=16)
pl.title(condname + ' Ear', fontsize=16)
pl.grid()
pl.xticks(fontsize=14)
pl.yticks(fontsize=14)


from scipy import signal
from anlffr.spectral import mtplv
params = dict(Fs=raw.info['sfreq'], tapers=[4, 7], fpass=[1, 4000], itc=0)
plv, f = mtplv(x[goods, :, :], params)
plv_mean = signal.savgol_filter(plv.mean(axis=0), 5, 3)

# Plot in SNR terms if needed
plotTMTF = False
pl.figure()
if plotTMTF:
    tmtf = plv_mean * ( (x.shape[1] * params['tapers'][1]))/ ( 1 - plv_mean)
    pl.plot(f, 10 * np.log10(tmtf), linewidth=2)
    pl.ylabel('TMTF SNR (dB)', fontsize=16)

else:
    pl.plot(f, plv_mean, linewidth=2)
    pl.ylabel('Phase Locking Value', fontsize=16)

pl.xlabel('Modulation Frequency (Hz)', fontsize=16)
pl.xticks(fontsize=14)
pl.yticks(fontsize=14)
pl.xscale('log')
pl.xticks([16, 32, 64, 128, 256, 512, 1024, 2048],
          labels=['16', '32', '64', '128', '256', '512', '1024', '2048'],
          fontsize=14)
pl.xlim((16, 2048))
pl.grid()


saveResults = True
if saveResults:
    mdict = dict(t=t, x=y, f=f, plv=plv)
    savepath = froot + '/CentralGainResults/'
    savename = subj + savestem + '_CentralGain.mat'
    savemat(savepath + savename, mdict)
    

