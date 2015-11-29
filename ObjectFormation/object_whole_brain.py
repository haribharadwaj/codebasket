"""
==================================================
Compute MNE-dSPM inverse solution on single epochs
==================================================

Compute dSPM inverse solution on single trial epochs restricted
to a brain label.

"""
# Author: Hari Bharadwaj <hari@nmr.mgh.harvard.edu>
#
# Compyright 2015. All Rights Reserved.

import numpy as np
import mne
from mne.minimum_norm import apply_inverse_epochs, read_inverse_operator
from anlffr.tfr import tfr_multitaper

froot = '/autofs/cluster/transcend/hari/ObjectFormation/'
subj = '093302'
para = 'object'
conds = ['onset', 'coh07']
sss = True
if sss:
    ssstag = '_sss'
else:
    ssstag = ''

froot_subj = froot + '/' + subj + '/'
fname_inv = froot_subj + subj + '_' + para + ssstag + '-inv.fif'

# Using the same inverse operator when inspecting single trials Vs. evoked
snr = 3.0  # Standard assumption for average data but using it for single trial
lambda2 = 1.0 / snr ** 2

method = "dSPM"  # use dSPM method (could also be MNE or sLORETA)

# Load operator and labels
inverse_operator = read_inverse_operator(fname_inv)

stc_results = []
for k, cond in enumerate(conds):

    # Read epochs
    fname_epochs = (froot_subj + subj + ssstag + '_' + para + '_' + cond +
                    '-epo.fif')
    epochs = mne.read_epochs(fname_epochs)
    times = epochs.times
    # Get evoked data (averaging across trials in sensor space)
    evoked = epochs.average()

    # Compute inverse solution and stcs for each epoch
    # Use the same inverse operator as with evoked data (i.e., set nave)
    # If you use a different nave, dSPM just scales by a factor sqrt(nave)
    stcs = apply_inverse_epochs(epochs, inverse_operator, lambda2, method,
                                label=None, nave=evoked.nave)

    #######################################################
    # Calculating whole-brain TFR stuff at only one frequency
    stc_result = stcs[0].copy()
    nVerts = stc_result.shape[0]
    for v in range(nVerts):
        print 'Doing vertex #', v,  '/', nVerts
        epo_array = np.zeros((len(stcs), 1, stcs[0].shape[1]))
        for k, stc in enumerate(stcs):
            epo_array[k, 0, :] = stc.data[v, :]

        freqs = [40., ]
        n_cycles = 9.
        power, itc, faketimes = tfr_multitaper(epo_array, epochs.info['sfreq'],
                                               freqs, n_cycles=n_cycles,
                                               time_bandwidth=2.0,
                                               zero_mean=True,
                                               verbose='DEBUG')
        stc_result.data[v, :] = power.squeeze()

    stc_stem = subj + '_' + para + '_gamma_' + cond
    stc_result.save(froot_subj + stc_stem)
    stc_results += [stc_result, ]

stc_contrast = stc_results[0] - stc_results[1]
stc_contrast_stem = subj + '_' + para + '_gamma_' + conds[0] + '-' + conds[1]
stc_contrast.save(froot_subj + stc_contrast_stem)
