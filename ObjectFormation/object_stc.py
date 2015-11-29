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
import pylab as pl
from mne.minimum_norm import apply_inverse_epochs, read_inverse_operator
from mne.minimum_norm import apply_inverse
from anlffr.tfr import tfr_multitaper, rescale, plot_tfr

froot = '/autofs/cluster/transcend/hari/ObjectFormation/'
subj = '093302'
para = 'object'
conds = ['coh20', 'coh14']
sss = True
if sss:
    ssstag = '_sss'
else:
    ssstag = ''

fname_inv = froot + '/' + subj + '/' + subj + '_' + para + ssstag + '-inv.fif'
label_name = 'object_manual-lh'
fname_label = froot + '/' + subj + '/' + subj + '_%s.label' % label_name

# Using the same inverse operator when inspecting single trials Vs. evoked
snr = 3.0  # Standard assumption for average data but using it for single trial
lambda2 = 1.0 / snr ** 2

method = "dSPM"  # use dSPM method (could also be MNE or sLORETA)

# Load operator and labels
inverse_operator = read_inverse_operator(fname_inv)
label = mne.read_label(fname_label)

powers = []
itcs = []

for k, cond in enumerate(conds):

    # Read epochs
    fname_epochs = (froot + '/' + subj + '/' + subj + ssstag + '_' + para +
                    '_' + cond + '-epo.fif')
    epochs = mne.read_epochs(fname_epochs)
    times = epochs.times
    # Get evoked data (averaging across trials in sensor space)
    evoked = epochs.average()

    # Compute inverse solution and stcs for each epoch
    # Use the same inverse operator as with evoked data (i.e., set nave)
    # If you use a different nave, dSPM just scales by a factor sqrt(nave)
    stcs = apply_inverse_epochs(epochs, inverse_operator, lambda2, method,
                                label, nave=evoked.nave)

    stc_evoked = apply_inverse(evoked, inverse_operator, lambda2, method)

    stc_evoked_label = stc_evoked.in_label(label)

    # Average over label (not caring to align polarities here)
    label_mean_evoked = np.mean(stc_evoked_label.data, axis=0)

    #######################################################
    # Calculating label TFR stuff

    epo_array = np.zeros((len(stcs), 1, stcs[0].shape[1]))
    nVerts = stc_evoked_label.shape[0]
    nTopVerts = 1
    corr_list = np.zeros(nVerts)
    for k in range(nVerts):
        corr_list[k] = np.corrcoef(label_mean_evoked,
                                   stc_evoked_label.data[k, :])[0, 1]
    topverts = np.argsort(corr_list)[-nTopVerts:]
    for k, stc in enumerate(stcs):
        epo_array[k, 0, :] = stc.data[topverts, :].mean(axis=0)

    freqs = np.arange(5., 130., 2.)
    n_cycles = freqs * 0.250
    power, itc, faketimes = tfr_multitaper(epo_array, epochs.info['sfreq'],
                                           freqs, n_cycles=n_cycles,
                                           time_bandwidth=2.0, zero_mean=False,
                                           verbose='DEBUG')
    powers += [power, ]
    itcs += [itcs, ]

power_contrast = 20. * np.log10(powers[0] / powers[1])
power_contrast_scaled = rescale(power_contrast, times, baseline=(-0.2, 0.),
                                mode='mean')
plot_tfr(power_contrast, times, freqs)
pl.show()
