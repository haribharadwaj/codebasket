from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch
import pylab as pl

# Adding Files and locations
# froot = '/cluster/transcend/hari/MEMR/'
froot = '/home/hari/Data/MEMR/'
# froot = '/Users/Hari/Documents/Data/MEMR/ABR/'
# froot = '/home/hari/Documents/PythonCodes/MEMR'

subjlist = ['I61', ]
ear = 'right'
MLR = False
if MLR:
    conds = [[3, 6]]
else:
    # conds = [[2, 5], [3, 6]]
    conds = [[3, 6]]

for subj in subjlist:

    fpath = froot + '/' + subj + '/'

    print 'Running Subject', subj
    rawlist = []
    evelist = []
    nruns = 6
    for run in range(nruns):
        bdfs = fnmatch.filter(os.listdir(fpath), subj + '_MEMR_' + ear +
                              '_' + str(run + 1) + '_ABR.bdf')

        if len(bdfs) == 1:
            edfname = fpath + bdfs[0]
            # Load data and read event channel
            if ear == 'right':
                refchans = ['EXG2']
            else:
                refchans = ['EXG1', ]
            (rawtemp, evestemp) = bs.importbdf(edfname, nchans=35,
                                               refchans=refchans)
            rawtemp.set_channel_types({'EXG3': 'eeg'})
            rawlist += [rawtemp, ]
            evelist += [evestemp, ]
        else:
            if len(bdfs) == 0:
                pass
            else:
                RuntimeError('More than one bdf file for same run!!')

    raw, eves = mne.concatenate_raws(rawlist, events_list=evelist)
    # Filter the data
    if MLR:
        raw.filter(l_freq=30., h_freq=1500, picks=np.arange(35))
        tmin, tmax = -0.005, 0.08
    else:
        raw.filter(l_freq=70., h_freq=3000, picks=np.arange(35))
        tmin, tmax = -0.002, 0.015
    raw.info['bads'] += ['EXG3', 'A1', 'A2', 'A30', 'A7', 'A6',
                         'A24', 'A28', 'A29', 'A3', 'A11', 'A15',
                         'A16', 'A17', 'A10', 'A21', 'A20', 'A25']
    goods = [28, 3, 30, 26, 4, 25, 7, 31, 22, 9, 8, 21, 11, 12, 18]
    abrs = []
    pl.figure()
    for cond in conds:
        print 'Doing condition ', cond
        epochs = mne.Epochs(raw, eves, cond, tmin=tmin, proj=False,
                            tmax=tmax, baseline=(0.001, 0.002),
                            picks=np.arange(35),
                            reject=dict(eeg=40e-6),
                            verbose='WARNING')
        abr = epochs.average()
        abrs += [abr, ]
        x = abr.data * 1e6  # microV
        t = abr.times * 1e3 - 1.0  # Adjust for delay and use milliseconds
        pl.plot(t, x[goods, :].mean(axis=0) - x[34, :], linewidth=2)
        pl.hold(True)
pl.xlabel('Time (ms)', fontsize=16)
if MLR:
    pl.ylabel('MLR (uV)', fontsize=16)
    pl.xlim((-5., 60.))
else:
    pl.ylabel('ABR (uV)', fontsize=16)
    pl.xlim((-4., 13.))

ax = pl.gca()
ax.tick_params(labelsize=16)
if not MLR:
    # pl.legend(('Condensation', 'Rarefaction'), loc='best')
    # pl.legend(('80 dB peSPL', '100 dB peSPL'), loc='best')
    pl.legend(('100 dB peSPL', ), loc='best')
pl.show()
mne.write_evokeds(fpath + subj + '_' + ear + '_strict-ave.fif', abrs)
