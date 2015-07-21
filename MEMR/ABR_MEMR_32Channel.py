from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch
import pylab as pl

# Adding Files and locations
# froot = '/cluster/transcend/hari/MEMR/'
# froot = '/Users/Hari/Documents/Data/MEMR/ABR/'
froot = '/home/hari/Documents/PythonCodes/MEMR'

subjlist = ['I02', ]
ear = 'Right'
conds = [[2, 5], [3, 6]]
for subj in subjlist:

    fpath = froot + '/' + subj + '/'

    print 'Running Subject', subj
    rawlist = []
    evelist = []
    nruns = 3
    for run in range(nruns):
        bdfs = fnmatch.filter(os.listdir(fpath), subj + '_MEMR_' + ear +
                              '_' + str(run + 1) + '_ABR.bdf')

        if len(bdfs) == 1:
            edfname = fpath + bdfs[0]
            # Load data and read event channel
            if ear == 'Right':
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
        # Filter the data
        # raw.filter(l_freq=33., h_freq=3000, picks=np.arange(30, 35, 1))
    raw, eves = mne.concatenate_raws(rawlist, events_list=evelist)
    abrs = []
    pl.figure()
    for cond in conds:
        print 'Doing condition ', cond
        epochs = mne.Epochs(raw, eves, cond, tmin=-0.005, proj=False,
                            tmax=0.014, baseline=(0., 0.002),
                            picks=np.arange(30, 35, 1),
                            reject=dict(eeg=100e-6),
                            verbose='WARNING')
        abr = epochs.average()
        abrs += [abr, ]
        x = abr.data * 1e6  # microV
        t = abr.times * 1e3 - 1.6  # Adjust for delay and use milliseconds
        pl.plot(t, x[1, :] - x[4, :], linewidth=2)
        pl.hold(True)
pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('ABR (uV)', fontsize=16)
pl.xlim((0., 9.))
ax = pl.gca()
ax.tick_params(labelsize=16)
# pl.legend(('Condensation', 'Rarefaction'), loc='best')
pl.legend(('80 dB peSPL', '100 dB peSPL'), loc='best')
pl.show()
mne.write_evokeds(fpath + subj + '_abr-ave.fif', abrs)
