# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 15:50:07 2018

@author: hari
"""

import numpy as np
import mne
import pylab as pl
from sklearn.decomposition import PCA
from scipy.signal import savgol_filter as sg


# subjlist = ['035201', '038301', '038302', '039001', '042201', '092002',
#            '096301', '096302', '053001', '030801', '032901', '032902',
#            '013703', '014002', '063101', '075401', '011302', '010401',
#            '096901', '096902', '097201', '097301', '097601', '097701',
#            '098001', '098002', '098101', '098501', '011201', '011202']

subjlist = []
nsubjs = len(subjlist)
coh07summary = np.zeros((nsubjs, 1401))
coh14summary = np.zeros((nsubjs, 1401))
coh20summary = np.zeros((nsubjs, 1401))

ncomps = 2
combinecomps = True

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

    start, stop = ref.time_as_index([-0.1, 0.5])
    x = ref.data[:, start:stop]
    tref = ref.times[start:stop]
    pca = PCA(n_components=ncomps)

    pca.fit(x.T)  # Estimate the model and store within object

    y = pca.transform(x.T)

    varexp = pca.explained_variance_ratio_.sum()
    print 'Variance Explained is', varexp * 100, '%'

    if varexp < 0.5:
        print 'Warning! Variance Explained is less than 60% !'

    if combinecomps:
        normfac = np.max((y ** 2.).sum(axis=1)) ** 0.5
    else:
        normfac = np.max(np.abs(y[:, 0]))

    t = ref.times
    coh07 = pca.transform(evokeds[0].data.T)
    coh14 = pca.transform(evokeds[1].data.T)
    coh20 = pca.transform(evokeds[2].data.T)

    if combinecomps:
        filtlen = 15
        filtord = 2
        coh07summary[k, :] = sg((coh07 ** 2.).sum(axis=1) ** 0.5 / normfac,
                                filtlen, filtord)
        coh14summary[k, :] = sg((coh14 ** 2.).sum(axis=1) ** 0.5 / normfac,
                                filtlen, filtord)
        coh20summary[k, :] = sg((coh20 ** 2.).sum(axis=1) ** 0.5 / normfac,
                                filtlen, filtord)
    else:
        coh07summary[k, :] = coh07[:, 0] / normfac
        coh14summary[k, :] = coh14[:, 0] / normfac
        coh20summary[k, :] = coh20[:, 0] / normfac

pl.figure()
m7 = coh07summary.mean(axis=0)
e7 = coh07summary.std(axis=0) / (nsubjs ** 0.5)
m14 = coh14summary.mean(axis=0)
e14 = coh14summary.std(axis=0) / (nsubjs ** 0.5)
m20 = coh20summary.mean(axis=0)
e20 = coh20summary.std(axis=0) / (nsubjs ** 0.5)
alpha = 0.25
pl.plot(t, m7)
pl.fill_between(t, m7 - e7, m7 + e7, alpha=alpha)
pl.plot(t, m14)
pl.fill_between(t, m14 - e14, m14 + e14, alpha=alpha)
pl.plot(t, m20)
pl.fill_between(t, m20 - e20, m20 + e20, alpha=alpha)
pl.legend(('7 / 20', '14 / 20', '20 / 20'), title='Number of Coherent Tones',
          fontsize=14)
pl.xlabel('Time (s)', fontsize=16)
pl.ylabel('Evoked Response (normalized)', fontsize=16)
pl.xlim((-0.15, 1.0))
ax = pl.gca()
ax.tick_params(labelsize=14)
pl.show()
