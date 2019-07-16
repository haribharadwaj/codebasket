from anlffr.helper import biosemi2mne as bs
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs
import mne
import os
import fnmatch
from mne.time_frequency import tfr_multitaper
import numpy as np
import pylab as pl

# Adding Files and locations
froot = '/home/hari/Data/Inhibition/'

subjlist = ['S160', ]

condlist = [[1, 5], [2, 6], [3, 7], [4, 8]]
condnames = ['nopre', 'same', 'low', 'high']
doITC = True
for subj in subjlist:
    evokeds = []
    itcs = []
    powers = []
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
                                         refchans=None, exclude=[])
        rawlist += [rawtemp, ]
        evelist += [evestemp, ]
    raw, eves = mne.concatenate_raws(rawlist, events_list=evelist)

    # raw.set_channel_types({'EXG3': 'eeg', 'EXG4': 'eeg'})
    raw.info['bads'] += []

    # Filter the data
    raw.filter(l_freq=1.5, h_freq=40.)

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
        tmin = -0.9
        epochs = mne.Epochs(
            raw, eves, cond, tmin=tmin, proj=True,
            tmax=1.8, baseline=(-0.5, 0.0),
            reject=dict(eeg=100e-6))
        evoked = epochs.average()
        evokeds += [evoked, ]

        if doITC:
            # Compute evoked response using ITC
            freqs = np.arange(1., 30., 1.)
            n_cycles = freqs * 0.2
            picks = [30, 31]
            power, itc = tfr_multitaper(epochs, freqs, n_cycles, picks=picks,
                                        time_bandwidth=4.0, n_jobs=-1)
            itc.apply_baseline(baseline=(-0.7, 0))
            power.apply_baseline(baseline=(-0.7, 0), mode='logratio')
            itcs += [itc, ]
            powers += [power, ]

    resname = fpath + subj + '_inh-ave.fif'
    mne.write_evokeds(resname, evokeds)

    if doITC:
        resnameitc = fpath + subj + '_inh_itc-tfr.h5'
        mne.time_frequency.write_tfrs(resnameitc, itcs, overwrite=True)
        resnamepow = fpath + subj + '_inh_pow-tfr.h5'
        mne.time_frequency.write_tfrs(resnamepow, powers, overwrite=True)

    # Plot single channel evoked responses for all conditions
    t = evoked.times
    x = np.zeros((t.shape[0], len(condnames)))
    ch = 31
    for k in range(len(condnames)):
        x[:, k] = evokeds[k].data[ch, :] * 1.0e6
    pl.plot(t, x)
    pl.xlabel('Time (s)')
    pl.ylabel('Evoked Response (uV)')
    pl.legend(condnames)

    t = itc.times  # Just in case
    fselect = freqs < 10.
    y = np.zeros((t.shape[0], len(condnames)))
    ch = 1
    for k in range(len(condnames)):
        y[:, k] = itcs[k].data[ch, fselect, :].squeeze().mean(axis=0)
    pl.figure()
    pl.plot(t, y)
    pl.xlabel('Time (s)')
    pl.ylabel('ITC (baseline subtracted)')
    pl.legend(condnames)
