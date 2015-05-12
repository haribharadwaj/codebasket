import mne
import numpy as np
import os
import fnmatch

from anlffr.helper import biosemi2mne as bs
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs
from mne.time_frequency import tfr_multitaper

# Adding Files and locations
# froot = '/Users/Hari/Documents/Data/ATTEEG/'
froot = '/Users/Hari/Documents/Data/ATTEEG/'

# List of files stems, each will be appended by run number
# Could be different for edf, fiff, eve etc.
# Use list [] and enumerate over if filenames are weird

subjlist = ['I33', ]

ch = [3, 4, 25, 26, 30, 31]  # Channels of interest for auditory ERPs

for k, subj in enumerate(subjlist):
    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    print 'Running Subject', subj

    bdfs = fnmatch.filter(os.listdir(fpath), subj + '*.bdf')

    if len(bdfs) > 1:
        print '***WARNING!! Multiple files found!***'

    print 'Viola! Data files found!'
    edfname = bdfs[0]

    # Load data and read event channel
    (raw, eves) = bs.importbdf(fpath + edfname, refchans=None)
    raw.info['bads'] += ['A7', 'A24', 'A28', 'EXG1', 'EXG2']

    # Filter the data for ERPs
    raw.filter(l_freq=0.5, h_freq=200, l_trans_bandwidth=0.15,
               picks=np.arange(0, 32, 1))
    # SSP for blinks
    blinks = find_blinks(raw)
    epochs_blinks = mne.Epochs(raw, blinks, 998, tmin=-0.25, tmax=0.25,
                               proj=True, baseline=(-0.25, 0),
                               reject=dict(eeg=500e-6))
    blink_projs = compute_proj_epochs(epochs_blinks, n_grad=0, n_mag=0,
                                      n_eeg=2, verbose='DEBUG')
    raw.add_proj(blink_projs)

    # Epoching and averaging
    epochsL = mne.Epochs(raw, eves, [3, ], tmin=-0.5, proj=True,
                         tmax=6.0, baseline=(-0.5, 0.0),
                         reject=dict(eeg=200e-6))
    epochsR = mne.Epochs(raw, eves, [5, ], tmin=-0.5, proj=True,
                         tmax=6.0, baseline=(-0.5, 0.0),
                         reject=dict(eeg=200e-6))

freqs = np.arange(2, 80, 2)  # define frequencies of interest
n_cycles = freqs / float(5)  # different number of cycle per frequency
n_cycles[freqs < 5] = 1
# epochs.subtract_evoked()
powerL, itcL = tfr_multitaper(epochsL, freqs=freqs, n_cycles=n_cycles,
                              time_bandwidth=2.0, n_jobs=-1)
powerL.data = 20*np.log10(powerL.data)
powerL.plot_topo(dB=False, baseline=(-0.5, 0.))

powerR, itcR = tfr_multitaper(epochsR, freqs=freqs, n_cycles=n_cycles,
                              time_bandwidth=2.0, n_jobs=-1)
powerR.data = 20*np.log10(powerR.data)
powerL.plot_topo(dB=False, baseline=(-0.5, 0.))

powerLR = powerR - powerL
powerLR.plot_topo(dB=False, baseline=(-0.5, 0.))
