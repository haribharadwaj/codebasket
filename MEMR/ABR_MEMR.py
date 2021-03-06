from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch
import pylab as pl

# Adding Files and locations
# froot = '/cluster/transcend/hari/MEMR/'
froot = '/Users/Hari/Documents/Data/MEMR/ABR/'

subjlist = ['I50', ]
ear = 'Right'
conds = [3, 6, [3, 6]]
for subj in subjlist:

    fpath = froot + '/' + subj + '/'

    print 'Running Subject', subj

    bdfs = fnmatch.filter(os.listdir(fpath), subj + '_' + ear +
                          '_ABR_module.bdf')

    if len(bdfs) > 1:
        print 'Warning! More than 1 file found!'
    else:
        edfname = bdfs[0]

    # Load data and read event channel
    (raw, eves) = bs.importbdf(fpath + edfname, nchans=2, refchans=None)

    # Filter the data
    raw.filter(l_freq=130, h_freq=2000, picks=np.arange(0, 2, 1))
    abrs = []
    pl.figure()
    for cond in conds:
        print 'Doing condition ', cond
        epochs = mne.Epochs(raw, eves, cond, tmin=0., proj=False,
                            tmax=0.014, baseline=(0.001, 0.002),
                            reject=dict(eeg=20e-6), verbose='WARNING')
        abr = epochs.average()
        abrs += [abr, ]
        x = abr.data * 1e6  # microV
        t = abr.times * 1e3 - 1.6  # Adjust for delay and use milliseconds
        if ear == 'Right':
            pl.plot(t, x[1, :] - x[0, :], linewidth=2)
        else:
            pl.plot(t, x[0, :] - x[1, :], linewidth=2)
        pl.hold(True)
pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('ABR (uV)', fontsize=16)
pl.legend(('Condensation', 'Rarefaction', 'Combined'), loc='best')
pl.xlim((0., 9.))
ax = pl.gca()
ax.tick_params(labelsize=16)
pl.show()
mne.write_evokeds(fpath + subj + '_abr-ave.fif', abrs)
