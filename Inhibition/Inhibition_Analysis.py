from anlffr.helper import biosemi2mne as bs
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs
import mne
import os
import fnmatch
from mne.time_frequency import tfr_multitaper
import numpy as np

# Adding Files and locations
froot = 'D:/DATA/Inhibition/'

subjlist = ['S160', ]

condlist = [[1, 5], [2, 6], [3, 7], [4, 8]]
condnames = ['nopre', 'same', 'low', 'high']
overwriteOld = True
for subj in subjlist:
    evokeds = []
    itcs = []
    # Load data and read event channel
    fpath = froot + subj + '/'
    bdfs = fnmatch.filter(os.listdir(fpath), subj +
                          '_inh*.bdf')
    print 'Viola!', len(bdfs),  'files found!'

    # Load data and read event channel
    rawlist = []
    evelist = []
    for k, rawname in enumerate(bdfs):
        rawtemp, evestemp = bs.importbdf(fpath + rawname, verbose='DEBUG',
                                         refchans=None)
        rawlist += [rawtemp, ]
        evelist += [evestemp, ]
    raw, eves = mne.concatenate_raws(rawlist, events_list=evelist)

    # raw.set_channel_types({'EXG3': 'eeg', 'EXG4': 'eeg'})
    raw.info['bads'] += []

    # Filter the data
    raw.filter(l_freq=1.5, h_freq=90.)

    removeblinks = True

    if removeblinks:
        # SSP for blinks
        blinks = find_blinks(raw, ch_name=['A1', ],  l_trans_bandwidth=0.4)
        blinkname = (fpath + subj + '_inh_blinks_erp' + '-eve.fif')
        mne.write_events(blinkname, blinks)
        epochs_blinks = mne.Epochs(raw, blinks, 998, tmin=-0.25,
                                   tmax=0.25, proj=True,
                                   baseline=(-0.25, 0),
                                   reject=dict(eeg=500e-6))
        blink_projs = compute_proj_epochs(epochs_blinks, n_grad=0,
                                          n_mag=0, n_eeg=2,
                                          verbose='DEBUG')
        raw.add_proj(blink_projs)

    for c, cond in enumerate(condlist):

        condname = condnames[c]
        # Epoching events of type
        epochs = mne.Epochs(
            raw, eves, cond, tmin=-0.3, proj=True,
            tmax=1.8, baseline=(-0.3, 0.0),
            reject=dict(eeg=200e-6))  # 200 regular, 50 strict
        evoked = epochs.average()
        evokeds += [evoked, ]
        freqs = np.arange(2., 80., 1.)
        n_cycles = freqs * 0.2
        power, itc = tfr_multitaper(epochs, freqs, n_cycles,
                                    time_bandwidth=2.0, n_jobs=-1)
        itc.apply_baseline(baseline=(-0.3, 0))
        itcs += [itc, ]

    resname = fpath + subj + '_inh-ave.fif'
    mne.write_evokeds(resname, evokeds)
    resnameitc = fpath + subj + '_inh_itc-tfr.h5'
    mne.time_frequency.write_tfrs(resnameitc, itcs)
