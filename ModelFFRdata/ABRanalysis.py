from anlffr.helper import biosemi2mne as bs
import mne
import numpy as np
import os
import fnmatch

# Adding Files and locations
froot = '/home/hari/Documents/PythonCodes/ModelData/'

# List of files stems, each will be appended by run number
# Could be different for edf, fiff, eve etc.
# Use list [] and enumerate over if filenames are weird


# subjlist =['I01','I02','I03','I06','I07','I08','I09','I11','I13','I14','I15',
# 'I17_redo','I18','I19','I20','I25','I26','I27','I28','I29','I30','I37',
# 'I16','I32','I33','I34','I35','I39','I05','I36']

subjlist = ['I14']

for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    cond = 1
    condstem = '_click'

    print 'Running Subject', subj, 'Condition', cond

    bdfs = fnmatch.filter(os.listdir(fpath), subj + '*ABR_10Hz.bdf')

    if len(bdfs) > 1:
        print 'Warning! More than 1 file found!'
    else:
        edfname = bdfs[0]

    # Load data and read event channel
    (raw, eves) = bs.importbdf(fpath + edfname, nchans=35,
                               refchans=['EXG1'])
    raw.info['bads'] += ['EXG3']
    # Filter the data
    raw.filter(
        l_freq=30, h_freq=2000, picks=np.arange(0, 34, 1))

    raw.notch_filter(freqs=np.arange(60, 2000, 60),
                     picks=np.arange(0, 34, 1))
    # Here click trigger is 1
    selectedEve = dict(up=cond)

    # Epoching events of type 1 and 7
    epochs = mne.Epochs(
        raw, eves, selectedEve, tmin=-0.02, proj=False,
        tmax=0.2, baseline=(-0.02, 0),
        reject = dict(eeg=50e-6))
    abr = epochs.average()
    abr.plot()
