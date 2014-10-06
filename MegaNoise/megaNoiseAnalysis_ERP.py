import mne
import numpy as np
import os
import fnmatch
from anlffr.helper import biosemi2mne as bs
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs
import pylab as pl
from mne.time_frequency.tfr import _induced_power as induced_power

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/MegaNoise/'

subjlist = ['I13', ]
nchans = 34

epochs = []

for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'
    switch = [True, False]
    evokeds = []
    for cond in switch:
        if cond:
            condlist = [11, 12, 21, 22, 33, 34, 43, 44]  # ITD trans
            # condlist = [11, 13, 31, 33, 22, 24, 42, 44]  # AM rate trans
            # condlist = [11, 13, 31, 33]  # 25-25
            # condlist = [11, 21, 31, 41, 12, 22, 32, 42]  # L
            # condlist = [11, 21, 12, 22]  # LL
            # condlist = [33, 34, 43, 44]  # RR
            condstem = 'Repeat'
        else:
            condlist = [13, 14, 23, 24, 31, 32, 41, 42]  # ITD trans
            # condlist = [12, 14, 21, 23, 32, 34, 41, 43]  # AM rate trans
            # condlist = [21, 41, 23, 43]  # 43-25
            # condlist = [13, 23, 33, 43, 14, 24, 34, 44]  # R
            # condlist = [31, 32, 41, 42]  # RL
            # condlist = [13, 14, 23, 24]  # LR
            condstem = 'Switch'

        print 'Running Subject', subj, 'Condition', condstem

        save_raw_name = subj + '_' + condstem + '_alltrial.mat'
        bdfs = fnmatch.filter(os.listdir(fpath), subj + '*.bdf')
        print 'Viola!', len(bdfs),  'files found!'
        if len(bdfs) > 1:
            print 'Wait!!.. Was expecting only one file..'
            'Going to use just one'
        edfname = bdfs[0]

        # Load data and read event channel
        (raw, eves) = bs.importbdf(fpath + edfname, nchans=nchans,
                                   refchans=['EXG1', 'EXG2'],
                                   verbose='DEBUG')

        raw.info['bads'] += ['EXG3', ]
        # Filter the data for ERPs
        raw.filter(l_freq=1.0, h_freq=20, l_trans_bandwidth=0.15,
                   picks=np.arange(0, 32, 1))

        # raw.apply_proj()
        fs = raw.info['sfreq']
        # SSP for blinks
        blinks = find_blinks(raw)
        epochs_blinks = mne.Epochs(raw, blinks, 998, tmin=-0.25,
                                   tmax=0.25, proj=True,
                                   baseline=(-0.25, 0),
                                   reject=dict(eeg=500e-6))
        blink_projs = compute_proj_epochs(epochs_blinks, n_grad=0,
                                          n_mag=0, n_eeg=2,
                                          verbose='DEBUG')
        # raw.del_proj(0)  # Removing average reference projection
        raw.add_proj(blink_projs)

        # Epoching events of type
        epochs = mne.Epochs(raw, eves, condlist, tmin=-0.3, proj=True,
                            tmax=1.0, baseline=(-0.3, 0.0),
                            reject = dict(eeg=100e-6))                   
        evokeds += [epochs.average(), ]

pl.figure()
ch = [3, 4, 25, 26, 30, 31]
ch = [4, ]
nconds = len(switch)
t = evokeds[0].times
x = evokeds[0].data[ch, :].mean(axis=0)*1e6
y = evokeds[1].data[ch, :].mean(axis=0)*1e6
z = x - y
pl.plot(t, x, linewidth=2)
pl.hold(True)
pl.plot(t, y, linewidth=2)
pl.plot(t, z, linewidth=2)
pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('[Switch - Repeat] Response (uV)', fontsize=16)
pl.legend(('Switch', 'Repeat', 'Difference'))
pl.axis('tight')
pl.show()

adapt = evokeds[0] - evokeds[1]
adapt.plot_topomap(times=[0.31, ], ch_type='eeg', vmin=-1, vmax=1)

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
