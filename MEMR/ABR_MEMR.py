from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch
import pylab as pl

# Adding Files and locations
froot = '/Users/Hari/Documents/Data/MEMR/'

subjlist = ['I13']
mne.set_log_level(verbose='WARNING')
for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    print 'Running Subject', subj

    bdfs = fnmatch.filter(os.listdir(fpath), subj + '*ABR*.bdf')

    if len(bdfs) > 1:
        print 'Warning! More than 1 file found!'
    else:
        edfname = bdfs[0]

    # Load data and read event channel
    (raw, eves) = bs.importbdf(fpath + edfname, nchans=2, refchans=None)

    # Filter the data
    raw.filter(l_freq=150, h_freq=2000, picks=np.arange(0, 2, 1))

    pre = range(1, 40)
    on = range(39, 79)
    post = range(79, 118)

    # Epoching events of type 3 and 8
    epochs_pre = mne.Epochs(raw, eves, pre, tmin=0., proj=False, tmax=0.014,
                            baseline=(0., 0.002), reject=dict(eeg=50e-6))
    abr_pre = epochs_pre.average()
    # Epoching events of type 3 and 8
    epochs_on = mne.Epochs(raw, eves, on, tmin=0., proj=False, tmax=0.014,
                           baseline=(0., 0.002), reject=dict(eeg=50e-6))
    abr_on = epochs_on.average()
    # Epoching events of type 3 and 8
    epochs_post = mne.Epochs(raw, eves, post, tmin=0., proj=False, tmax=0.014,
                             baseline=(0., 0.002), reject=dict(eeg=50e-6))
    abr_post = epochs_post.average()

x_pre = abr_pre.data * 1e6  # microV
t = abr_pre.times * 1e3 - 1.6  # Adjust for delay and use milliseconds
x_on = abr_on.data * 1e6  # microV
x_post = abr_post.data * 1e6  # microV
pl.figure()
pl.plot(t, x_pre[0, :] - x_pre[1, :], linewidth=2)
pl.hold(True)
pl.plot(t, x_on[0, :] - x_on[1, :], linewidth=2)
pl.hold(True)
pl.plot(t, x_post[0, :] - x_post[1, :], linewidth=2)
pl.xlabel('Time (ms)')
pl.ylabel('ABR (uV)')
pl.legend(('PRE', 'NOISE', 'POST'), loc=0)
pl.show()
