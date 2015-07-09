import numpy as np
import os
import fnmatch
import pylab as pl
from scipy.io import loadmat
from anlffr.preproc import band_pass_filter, peak_finder
from anlffr.spectral import mtspec
# Adding Files and locations
froot = '/Users/Hari/Documents/Data/MEMR/'

subjlist = ['I13']
fs = 48828.125  # Hz
OAE = True
for subj in subjlist:

    fpath = froot + subj + '/OAE/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    print 'Running Subject', subj

    matnames = fnmatch.filter(os.listdir(fpath), subj + '*CEOAE_MEMR*.mat')
    if not OAE:
        for kfile, matfile in enumerate(matnames):
            Praw = loadmat(fpath + matfile)['Pcanal'].squeeze()
            P = band_pass_filter(Praw, fs, 250, 20e3, filter_length='5ms')
            locs, peaks = peak_finder(P, thresh=0.2)
            if locs.shape[0] < 117 or locs.shape[0] > 115:
                nclicks_current = locs.shape[0]
                click_noise = np.zeros((39, 150))
                for k, loc in enumerate(locs):
                    if k < 78 and k > 38:
                        ind = np.arange(-50, 100) + loc
                        click_noise[k - 39, :] = P[ind]
                click_pre = np.zeros((39, 150))
                for k, loc in enumerate(locs):
                    if k < 39:
                        ind = np.arange(-50, 100) + loc
                        click_pre[k, :] = P[ind]
            if kfile == 0:
                pre = click_pre
                noise = click_noise
            else:
                pre = np.concatenate((pre, click_pre), axis=0)
                noise = np.concatenate((noise, click_noise), axis=0)
    else:
        for kfile, matfile in enumerate(matnames):
            Praw = loadmat(fpath + matfile)['Pcanal'].squeeze()
            P = band_pass_filter(Praw, fs, 250, 20e3, filter_length='5ms')
            locs, peaks = peak_finder(P, thresh=0.2)
            if locs.shape[0] < 117 or locs.shape[0] > 115:
                nclicks_current = locs.shape[0]
                click_noise = np.zeros((39, 900))
                for k, loc in enumerate(locs):
                    if k < 78 and k > 38:
                        ind = np.arange(100, 1000) + loc
                        click_noise[k - 39, :] = P[ind]
                click_pre = np.zeros((39, 900))
                for k, loc in enumerate(locs):
                    if k < 39:
                        ind = np.arange(100, 1000) + loc
                        click_pre[k, :] = P[ind]
            if kfile == 0:
                pre = click_pre
                noise = click_noise
            else:
                pre = np.concatenate((pre, click_pre), axis=0)
                noise = np.concatenate((noise, click_noise), axis=0)

params = dict(Fs=fs, tapers=[2, 3], fpass=[250., 20500.])
S, N, f = mtspec(pre - noise, params)
f = f / 1e3
pl.plot(f, np.log10(S) * 20, 'b', linewidth=2)
pl.hold(True)
pl.semilogx(f, np.log10(N) * 20, 'r--', linewidth=2)
pl.xlim((0.5, 20.5))
pl.xlabel('Frequency (kHz)', fontsize=16)
pl.ylabel('CEOAE (dB)', fontsize=20)
