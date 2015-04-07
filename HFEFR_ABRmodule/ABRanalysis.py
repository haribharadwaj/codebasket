from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch
import pylab as pl

# Adding Files and locations
froot = '/Users/Hari/Documents/Data/ABRmoduleBiosemi/'

subjlist = ['I13']

for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    condlist = [3, 8]
    condstemlist = ['_clickpos', '_clickneg']

    print 'Running Subject', subj

    bdfs = fnmatch.filter(os.listdir(fpath), subj + '*ABR*.bdf')

    if len(bdfs) > 1:
        print 'Warning! More than 1 file found!'
    else:
        edfname = bdfs[0]

    # Load data and read event channel
    (raw, eves) = bs.importbdf(fpath + edfname, nchans=2)

    # Filter the data
    raw.filter(l_freq=100, h_freq=3000, picks=np.arange(0, 2, 1))

    selectedEve = dict(pos=condlist[0], neg=condlist[1])

    # Epoching events of type 3 and 8
    epochs = mne.Epochs(
        raw, eves, selectedEve, tmin=0., proj=False,
        tmax=0.014, baseline=(0., 0.002), reject=dict(eeg=20e-6))
    abr = epochs.average()

x = abr.data * 1e6  # microV
t = abr.times * 1e3 - 1.6  # Adjust for delay and use milliseconds
pl.figure()
pl.plot(t, x[0, :] - x[1, :], linewidth=2)
pl.hold(True)
pl.plot(t, x[0, :], linewidth=2)
pl.xlabel('Time (ms)')
pl.ylabel('ABR (uV)')
pl.legend(('ipsi - contra ear', 'ipsi ear - midline'), loc=0)
pl.show()
