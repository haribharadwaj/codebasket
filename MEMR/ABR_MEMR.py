from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch
import pylab as pl

# Adding Files and locations
froot = '/cluster/transcend/hari/MEMR/'

subjlist = ['I14']
conds = [1, 2, 3, 4, 5, 6]
for subj in subjlist:

    fpath = froot + '/'

    print 'Running Subject', subj

    bdfs = fnmatch.filter(os.listdir(fpath), subj + '*ABR*.bdf')

    if len(bdfs) > 1:
        print 'Warning! More than 1 file found!'
    else:
        edfname = bdfs[0]

    # Load data and read event channel
    (raw, eves) = bs.importbdf(fpath + edfname, nchans=2, refchans=None)

    # Filter the data
    raw.filter(l_freq=100, h_freq=3000, picks=np.arange(0, 2, 1))
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
        pl.plot(t, x[0, :] - x[1, :], linewidth=2)
        pl.hold(True)
pl.xlabel('Time (ms)', fontsize=16)
pl.ylabel('ABR (uV)', fontsize=16)
pl.legend(('70 dB peSPL', '85 dB peSPL', '100 dB peSPL'), loc='best')
pl.xlim((0., 9.))
ax = pl.gca()
ax.tick_params(labelsize=16)
pl.show()
mne.write_evokeds(fpath + subj + '_abr-ave.fif', abrs)
