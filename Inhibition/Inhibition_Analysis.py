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

subjlist = ['S001', ]

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

    if subj == 'S001':
        exclude = [u'B1', u'B2', u'B3', u'B4', u'B5', u'B6', u'B7', u'B8',
                   u'B9', u'B10', u'B11', u'B12', u'B13', u'B14', u'B15',
                   u'B16', u'B17', u'B18', u'B19', u'B20', u'B21', u'B22',
                   u'B23', u'B24', u'B25', u'B26', u'B27', u'B28', u'B29',
                   u'B30', u'B31', u'B32', u'C1', u'C2', u'C3', u'C4', u'C5',
                   u'C6', u'C7', u'C8', u'C9', u'C10', u'C11', u'C12', u'C13',
                   u'C14', u'C15', u'C16', u'C17', u'C18', u'C19', u'C20',
                   u'C21', u'C22', u'C23', u'C24', u'C25', u'C26', u'C27',
                   u'C28', u'C29', u'C30', u'C31', u'C32', u'D1', u'D2',
                   u'D3', u'D4', u'D5', u'D6', u'D7', u'D8', u'D9', u'D10',
                   u'D11', u'D12', u'D13', u'D14', u'D15', u'D16', u'D17',
                   u'D18', u'D19', u'D20', u'D21', u'D22', u'D23', u'D24',
                   u'D25', u'D26', u'D27', u'D28', u'D29', u'D30', u'D31',
                   u'D32']
    else:
        exclude = []

    for k, rawname in enumerate(bdfs):
        rawtemp, evestemp = bs.importbdf(fpath + rawname, verbose='DEBUG',
                                         refchans=['EXG1', 'EXG2'],
                                         exclude=exclude)
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
            itc.apply_baseline(baseline=(-0.5, 0))
            power.apply_baseline(baseline=(-0.5, 0), mode='logratio')
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
    ch = 30
    for k in range(len(condnames)):
        x[:, k] = evokeds[k].data[ch, :] * 1.0e6
    pl.plot(t, x)
    pl.xlabel('Time (s)')
    pl.ylabel('Evoked Response (uV)')
    pl.legend(condnames)

    t = itc.times  # Just in case
    fselect = freqs < 15.
    y = np.zeros((t.shape[0], len(condnames)))
    for k in range(len(condnames)):
        perChan = itcs[k].data[:, fselect, :].mean(axis=1)
        y[:, k] = perChan.mean(axis=0) ** 2.
    pl.figure()
    pl.plot(t, y)
    pl.xlabel('Time (s)')
    pl.ylabel('ITC (baseline subtracted)')
    pl.xlim((-0.5, 1.5))
    pl.legend(condnames)
