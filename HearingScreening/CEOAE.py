import numpy as np
import os
import fnmatch
import pylab as pl
from scipy.io import loadmat
from anlffr.preproc import band_pass_filter, peak_finder
from anlffr.dpss import dpss_windows
# Adding Files and locations
froot = '/autofs/cluster/transcend/hari/HearingScreening/CEOAE/Results/'

subjlist = ['Manny_Right']
fs = 48828.125  # Hz

for subj in subjlist:

    fpath = froot + subj + '/'

    print 'Running Subject', subj

    matnames = fnmatch.filter(os.listdir(fpath), subj + '*CEOAE*.mat')
    for kfile, matfile in enumerate(matnames):
        inputClick = loadmat(fpath + matfile)['y'].squeeze()
        Praw = loadmat(fpath + matfile)['Pcanal'].squeeze()
        P = band_pass_filter(Praw, fs, 250, 20e3, filter_length='5ms')
        locs, peaks = peak_finder(inputClick, thresh=0.2)
        if locs.shape[0] < 117 or locs.shape[0] > 115:
            nclicks_current = locs.shape[0]
            click = np.zeros((nclicks_current, 900))
            for k, loc in enumerate(locs):
                ind = np.arange(250, 1150) + loc
                click[k, :] = P[ind]
        if kfile == 0:
            clicks = click
        else:
            clicks = np.concatenate((clicks, click), axis=0)

ceoae = np.median(clicks, axis=0).squeeze()
clicks_noise = clicks
clicks_noise[::2, ] *= -1.0
noise = np.median(clicks_noise, axis=0).squeeze()
f = np.arange(ceoae.shape[0]) * fs / ceoae.shape[0]
w, e = dpss_windows(ceoae.shape[0], 1., 1.)
S = np.fft.fft(w * ceoae).squeeze()
N = np.fft.fft(w * noise).squeeze()
f = f / 1e3
pl.plot(f, np.log10(np.abs(S)) * 20, 'b', linewidth=2)
pl.hold(True)
pl.plot(f, np.log10(np.abs(N)) * 20, 'r--', linewidth=2)
pl.xlim((0.5, 8.0))
pl.xlabel('Frequency (kHz)', fontsize=16)
pl.ylabel('CEOAE (dB)', fontsize=20)
pl.show()
