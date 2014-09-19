import mne
import numpy as np
import os
import fnmatch
from anlffr.helper import biosemi2mne as bs
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs
import pylab as pl

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/MegaNoise/'

lowLevel = False
if lowLevel:
    froot = froot + 'Level58dB/'

subjlist = ['I13', ]
nchans = 34
ch = [3, 4, 25, 26, 30, 31]  # Channels of interest
freqs = np.arange(5, 500, 2)  # define frequencies of interest
n_cycles = freqs / float(3)  # different number of cycle per frequency
n_cycles[freqs < 15] = 2

for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'
    switch = [False, True]
    evokeds = []
    for cond in switch:
        if cond:
            condlist = [11, 12, 21, 22, 33, 34, 43, 44]
            condstem = 'repeat'
        else:
            condlist = [13, 14, 23, 24, 31, 32, 41, 42]
            condstem = 'switch'

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
        raw.filter(l_freq=0.2, h_freq=50, l_trans_bandwidth=0.15,
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
        epochs = mne.Epochs(raw, eves, condlist, tmin=-0.2, proj=True,
                            tmax=0.45, baseline=(-0.02, 0.02),
                            reject = dict(eeg=150e-6))
        evokeds += [epochs.average(), ]

adapt = evokeds[1] - evokeds[0]
pl.figure()
ch = [3, 4, 25, 26, 30, 31]
nconds = len(switch)
t = evokeds[0].times
x = adapt.data[ch, :].mean(axis=0)*1e6
pl.plot(t, x, linewidth=2)
pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('[Switch - Repeat] Response (uV)', fontsize=16)
pl.show()
adapt.plot_topomap(times=[0.165, ], ch_type='eeg', vmin=-1, vmax=1)
