# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 15:50:07 2018

@author: hari
"""

import numpy as np
import mne
import pylab as pl
from sklearn.decomposition import PCA
from scipy import io


def mad(data, axis=None):
    return np.median(np.abs(data - np.median(data, axis)), axis)

tdlist = ['011201', '011202', '011302', '013703', '014002', '032901',
          '032902', '038301', '038302', '039001', '042201', '052402',
          '052901', '052902', '082601', '082802', '089401', '089402',
          '092002', '093101', '093302', '096301', '096302', '096603',
          '096901', '096902']


asdlist = ['010401', '030801', '035201', '053001', '063101', '075401',
           '082901', '085701', '086901', '087401', '092301', '093901',
           '095801', '097201', '097301', '097601', '097701', '098001',
           '098002', '098101', '098501']

tdage = [16, 15, 15, 15, 16, 14, 11, 12, 10, 15, 17, 10, 10, 9, 11, 7, 17,
         15, 16, 9, 17, 17, 8, 10, 15, 12]

asdage = [16, 17, 12, 12, 15, 12, 10, 15, 17, 17, 14, 10, 7, 16, 12, 15, 16,
          12, 12, 13, 16]

subjlist = tdlist + asdlist

age = tdage + asdage

nsubjs = len(subjlist)
coh07summary = np.zeros((nsubjs, 65))
coh14summary = np.zeros((nsubjs, 65))
coh20summary = np.zeros((nsubjs, 65))
normfacs = np.zeros(nsubjs)
varexps = np.zeros(nsubjs)
ncomps = 4
saveRes = True
zscore = True
for k, subj in enumerate(subjlist):
    print 'Loading subject', subj
    fname = './' + subj + '_sss_object_collapse-ave.fif'
    evokeds = mne.read_evokeds(fname, verbose='WARNING')
    ref = evokeds[3]  # Use the collapsed onset as reference

    if ref.info['sfreq'] == 3000:
        decim = 3
        ref.decimate(decim)
        evokeds[0].decimate(decim)
        evokeds[1].decimate(decim)
        evokeds[2].decimate(decim)

    if zscore:
        bstart, bstop = ref.time_as_index([-0.25, 0.])
        bmean = ref.data[:, bstart:bstop].mean(axis=1)
        bstd = ref.data[:, bstart:bstop].std(axis=1)
        ref.data = (ref.data.T - bmean).T
        ref.data = (ref.data.T / bstd).T

    start, stop = ref.time_as_index([-0.4, 1.0])

    x = ref.data[:, start:stop]
    tref = ref.times[start:stop]
    pca = PCA(n_components=ncomps)

    pca.fit(x.T)  # Estimate the model and store within object

    wts = (pca.components_ ** 2).sum(axis=0)
    wts[wts > max(wts)/2] = 0

    # Read TFR power data for 20 coherence only
    tfrname = './' + subj + '_sss_object_pow_coh20-tfr.h5'
    power = mne.time_frequency.read_tfrs(tfrname)
    x = power[0].data
    t = power[0].times
    x = x.transpose((1, 2, 0))
    powsummary = 10*np.log10((x[:, t > 0., :] * wts).sum(axis=2).mean(axis=1))
    powsummary -= max(powsummary)
    coh20summary[k, :] = powsummary

    # Read TFR power data for 20 coherence only
    tfrname = './' + subj + '_sss_object_pow_coh14-tfr.h5'
    power = mne.time_frequency.read_tfrs(tfrname)
    x = power[0].data
    t = power[0].times
    x = x.transpose((1, 2, 0))
    powsummary = 10*np.log10((x[:, t > 0., :] * wts).sum(axis=2).mean(axis=1))
    powsummary -= max(powsummary)
    coh14summary[k, :] = powsummary

    # Read TFR power data for 20 coherence only
    tfrname = './' + subj + '_sss_object_pow_coh07-tfr.h5'
    power = mne.time_frequency.read_tfrs(tfrname)
    x = power[0].data
    t = power[0].times
    x = x.transpose((1, 2, 0))
    powsummary = 10*np.log10((x[:, t > 0., :] * wts).sum(axis=2).mean(axis=1))
    powsummary -= max(powsummary)
    coh07summary[k, :] = powsummary
    f = power[0].freqs

m7 = coh07summary.mean(axis=0)
e7 = coh07summary.std(axis=0) / (nsubjs ** 0.5)
m14 = coh14summary.mean(axis=0)
e14 = coh14summary.std(axis=0) / (nsubjs ** 0.5)
m20 = coh20summary.mean(axis=0)
e20 = coh20summary.std(axis=0) / (nsubjs ** 0.5)
alpha = 0.25
pl.plot(f, m7)
pl.fill_between(f, m7 - e7, m7 + e7, alpha=alpha)
pl.plot(f, m14)
pl.fill_between(f, m14 - e14, m14 + e14, alpha=alpha)
pl.plot(f, m20)
pl.fill_between(f, m20 - e20, m20 + e20, alpha=alpha)
pl.legend(('7 / 20', '14 / 20', '20 / 20'),
          title='Number of Coherent Tones',
          fontsize=14)
pl.xlabel('Time (s)', fontsize=16)
pl.ylabel('Spectral Power (normalized)', fontsize=16)
ax = pl.gca()
ax.tick_params(labelsize=14)
pl.show()

if saveRes:
    mdict = dict(c7=coh07summary, c14=coh14summary, c20=coh20summary, f=f,
                 subjlist=subjlist, ntd=26, nasd=21, age=age)
    io.savemat('TFRPowerSummary.mat', mdict)
