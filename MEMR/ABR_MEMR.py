from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch
import pylab as pl

# Adding Files and locations
froot = '/Users/Hari/Documents/Data/MEMR/'

subjlist = ['I14']
conds = [4, 5]
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
    raw.filter(l_freq=100, h_freq=2000, picks=np.arange(0, 2, 1))
    abrs = []
    pl.figure()
    for cond in conds:
        epochs = mne.Epochs(raw, eves, cond, tmin=-0.001, proj=False,
                            tmax=0.014, baseline=(-0.001, 0.001),
                            reject=dict(eeg=25e-6), verbose='INFO')
        abr = epochs.average()
        abrs += [abr, ]
        x = abr.data * 1e6  # microV
        t = abr.times * 1e3 - 1.6  # Adjust for delay and use milliseconds
        pl.plot(t, x[1, :] - x[0, :], linewidth=2)
        pl.hold(True)
pl.xlabel('Time (ms)')
pl.ylabel('ABR (uV)')
# pl.legend(('60 dB peSPL', '70 dB peSPL', '80 dB peSPL', '90 dB peSPL',
#          '100 dB peSPL'))
pl.legend(('90 dB peSPL', '100 dB peSPL'))
pl.show()
mne.write_evokeds(fpath + subj + '_abr_conds4_5-ave.fif', abrs)
