from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch
import pylab as pl

# Adding Files and locations
froot = 'D:/DATA/ABR/'

subjlist = ['S012', ]

conds = [[3, 9], [5, 10], [6, 12], [48, 144], [80, 160], [96, 192]]

for subj in subjlist:

    fpath = froot + '/' + subj + '/'

    print 'Running Subject', subj
    rawlist = []
    evelist = []

    bdfs = fnmatch.filter(os.listdir(fpath), subj + '_ABR*.bdf')

    if len(bdfs) >= 1:
        for k, bdf in enumerate(bdfs):
            edfname = fpath + bdf
            # Load data and read event channel
            (rawtemp, evestemp) = bs.importbdf(edfname, nchans=36)
            rawtemp.set_channel_types({'EXG3': 'eeg', 'EXG4': 'eeg'})
            rawlist += [rawtemp, ]
            evelist += [evestemp, ]
    else:
        RuntimeError('No BDF files found!!')

    raw, eves = mne.concatenate_raws(rawlist, events_list=evelist)
    # Filter the data
    raw.filter(l_freq=70., h_freq=3000, picks=np.arange(36))
    tmin, tmax = -0.002, 0.015
    raw.info['bads'] += ['EXG3', 'EXG4', 'A1', 'A2', 'A30', 'A7', 'A6',
                         'A24', 'A28', 'A29', 'A3', 'A11', 'A15',
                         'A16', 'A17', 'A10', 'A21', 'A20', 'A25']
    goods = [28, 3, 30, 26, 4, 25, 7, 31, 22, 9, 8, 21, 11, 12, 18]
    abrs = []
    pl.figure()
    for cond in conds:
        print 'Doing condition ', cond
        epochs = mne.Epochs(raw, eves, cond, tmin=tmin, proj=False,
                            tmax=tmax, baseline=(0.001, 0.002),
                            picks=np.arange(36),
                            reject=dict(eeg=70e-6),
                            verbose='WARNING')
        abr = epochs.average()
        abrs += [abr, ]

x = abr.data * 1e6  # microV
t = abr.times * 1e3 - 1.0  # Adjust for delay and use milliseconds
pl.plot(t, x[goods, :].mean(axis=0) - x[35, :], linewidth=2)
pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('ABR (uV)', fontsize=16)
pl.xlim((0.5, 10.))
ax = pl.gca()
ax.tick_params(labelsize=16)

pl.show()
mne.write_evokeds(fpath + subj + '_ABR-ave.fif', abrs)
