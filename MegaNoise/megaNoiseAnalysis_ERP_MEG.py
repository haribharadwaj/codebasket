import mne
import numpy as np
import os
import fnmatch
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs
import pylab as pl
from mne.time_frequency.tfr import _induced_power as induced_power

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/ASSRmartinos/'

subjlist = ['SK', ]

epochs = []

for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'
    switch = [True, False]

    fifs = fnmatch.filter(os.listdir(fpath), subj + '*raw.fif')
    print 'Viola!', len(fifs),  'files found!'
    if len(fifs) > 1:
        print 'Wait!!.. Was expecting only one file..'
        'Going to use just one'
    fifname = fifs[0]

    # Load data and read event channel
    raw = mne.io.Raw(fpath + fifname, preload=True)
    eves = mne.find_events(raw, stim_channel='STI101', shortest_event=1)

    raw.info['bads'] += ['MEG1421', 'MEG1431', 'MEG2621']
    # Filter the data for ERPs
    raw.filter(l_freq=1.0, h_freq=20, l_trans_bandwidth=0.15,
               picks=np.arange(0, 308, 1))

    # raw.apply_proj()
    fs = raw.info['sfreq']
    # SSP for blinks
    blinks = find_blinks(raw, ch_name='EOG061')
    epochs_blinks = mne.Epochs(raw, blinks, 998, tmin=-0.25,
                               tmax=0.25, proj=True,
                               baseline=(-0.25, 0),
                               reject=dict(grad=8000e-13,
                                           mag=8e-12))
    blink_projs = compute_proj_epochs(epochs_blinks, n_grad=2,
                                      n_mag=2, n_eeg=0,
                                      verbose='DEBUG')
    raw.add_proj(blink_projs)
    evokeds = []

    for cond in switch:
        if cond:
            condlist = [1, 2, 5, 6, 11, 12, 15, 16]  # ITD trans
            condstem = 'Repeat'
        else:
            condlist = [3, 4, 7, 8, 9, 10, 13, 14]  # ITD trans
            condstem = 'Switch'

        print 'Running Subject', subj, 'Condition', condstem

        # Epoching events of type
        epochs = mne.Epochs(raw, eves, condlist, tmin=-0.3, proj=True,
                            tmax=1.0, baseline=(-0.3, 0.0),
                            reject=dict(grad=5000e-13, mag=5e-12))
        evokeds += [epochs.average(), ]

pl.figure()
ch = [175, ]
nconds = len(switch)
t = evokeds[0].times
x = evokeds[0].data[ch, :].mean(axis=0)*1e12
y = evokeds[1].data[ch, :].mean(axis=0)*1e12
z = x - y
pl.plot(t, x, linewidth=2)
pl.hold(True)
pl.plot(t, y, linewidth=2)
pl.plot(t, z, linewidth=2)
pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('[Switch - Repeat] Response (pT/cm)', fontsize=16)
pl.legend(('Switch', 'Repeat', 'Difference'))
pl.axis('tight')
pl.show()

adapt = evokeds[0] - evokeds[1]
adapt.plot_topomap(times=[0.1, 0.21], ch_type='grad')

##############################################################################
# ITC Analysis

freqs = np.arange(2, 100, 2)  # define frequencies of interest
n_cycles = freqs / float(5)  # different number of cycle per frequency
n_cycles[freqs < 5] = 1
dat1 = evokeds[0][:, ch, :].mean(axis=1, keepdims=True)
pow1, plv1 = induced_power(dat1, sfreq=fs, frequencies=freqs,
                           n_cycles=n_cycles, zero_mean=True)
dat2 = evokeds[1][:, ch, :].mean(axis=1, keepdims=True)
pow2, plv2 = induced_power(dat2, sfreq=fs, frequencies=freqs,
                           n_cycles=n_cycles, zero_mean=True)
fmin = 4
fmax = 10
bmin = -0.3
bmax = 0.0
plv1 = plv1.squeeze()
plv1 = plv1[np.logical_and(freqs > fmin, freqs < fmax), :].mean(axis=0)
plv2 = plv2[:, np.logical_and(freqs > fmin, freqs < fmax), :].mean(axis=2)
