import mne
import numpy as np
import os
import fnmatch
from anlffr.preproc import find_blinks


def pow2db(x):
    """ Converts *power* to decibels

    Parameters
    ----------

    x - Input in linear units

    Returns
    -------

    y - Equivalend in decibel units
    """

    y = 10*np.log10(x)
    return y

# Adding Files and locations
# froot = '/autofs/cluster/transcend/hari/MEMR/'
froot = '/Users/Hari/Documents/Data/ASSR_MEG/'
saveResults = True
subjlist = ['HB', ]

ch = range(1, 307)  # Channels of interest
mags = range(2, 306, 3)
grads = range(0, 306, 3) + range(1, 306, 3)

for subj in subjlist:
    fpath = froot + subj + '/'

    print 'Running Subject', subj

    save_raw_name = subj + '_ABR-epo.fif'

    if os.path.isfile(fpath + save_raw_name):
        preEpoched = True
        print 'Epoched data is already available on disk!'
        print 'Loading data from:', fpath + save_raw_name
        epochs = mne.read_epochs(fpath + save_raw_name, verbose='DEBUG')
        Fs = epochs.info['sfreq']
        x = epochs.get_data()
        times = epochs.times
    else:
        preEpoched = False
        fifs = fnmatch.filter(os.listdir(fpath), subj + '_CABR_raw.fif')
        print 'No pre-epoched data found, looking for raw files'
        print 'Viola!', len(fifs),  'files found!'
        for k, fif in enumerate(fifs):
            fifs[k] = fpath + fif
        # Load data and read event channel
        raw = mne.io.Raw(fifs, preload=True, add_eeg_ref=False)
        raw.set_channel_types({'EMG061': 'eeg'})
        raw.info['bads'] += ['MEG0223', 'MEG1623']

        # Filter the data for SSRs
        raw.filter(l_freq=70., h_freq=None, picks=None)
        eves = find_blinks(raw, event_id=1, thresh=0.1, l_freq=70.,
                           h_freq=2499., ch_name=['MISC001'])
        cond = 1

        # Epoching events of type
        epochs = mne.Epochs(raw, eves, cond, tmin=-0.005, proj=False,
                            tmax=0.050, baseline=(-0.005, 0.005),
                            reject=dict(grad=1000e-12, mag=10e-12),
                            verbose='WARNING')
abr = epochs.average()
avename = subj + '_ABR-ave.fif'
abr.save(fpath + avename)
