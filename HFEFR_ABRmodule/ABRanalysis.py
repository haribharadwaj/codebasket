from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch
import pylab as pl

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/HFEFR/'

# List of files stems, each will be appended by run number
# Could be different for edf, fiff, eve etc.
# Use list [] and enumerate over if filenames are weird


# subjlist =['I01','I02','I03','I06','I07','I08','I09','I11','I13','I14','I15',
# 'I17_redo','I18','I19','I20','I25','I26','I27','I28','I29','I30','I37',
# 'I16','I32','I33','I34','I35','I39','I05','I36']

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
    raw.filter(
        l_freq=200, h_freq=3000, picks=np.arange(0, 2, 1))

    #raw.notch_filter(freqs=np.arange(60, 3000, 60),
    #                 picks=np.arange(0, 2, 1))
    # Here click trigger is 1
    selectedEve = dict(pos=condlist[0], neg=condlist[1])

    # Epoching events of type 3 and 8
    epochs = mne.Epochs(
        raw, eves, selectedEve, tmin=-0.002, proj=False,
        tmax=0.014, baseline=(-0.002, 0),
        reject = dict(eeg=20e-6))
    abr = epochs.average()
    # abr.plot()

x = abr.data * 1e7
t = abr.times * 1e3 - 1.6  # Adjust for delay
pl.figure()
pl.plot(t, x[1, :], linewidth=2)
pl.hold(True)
pl.plot(t, x[0, :], linewidth=2)
pl.xlabel('Time (ms)')
pl.ylabel('ABR (uV)')
pl.legend(('Midline - contra ear', 'Midline - ipsi ear'), loc=0)
pl.show()
