import mne
import numpy as np
import os
import fnmatch
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs


# Adding Files and locations
# froot = '/home/hari/Documents/PythonCodes/ASSRmartinos/'
froot = '/autofs/cluster/transcend/hari/ASSRnew/'

subjlist = ['082601', ]

epochs = []

for subj in subjlist:

    fpath = froot + subj + '/'

    # These are so that the generated files are organized better
    respath = fpath + 'RES/'

    fifs = fnmatch.filter(os.listdir(fpath), subj + '*raw.fif')
    print 'Viola!', len(fifs),  'files found!'
    if len(fifs) > 1:
        print 'Wait!!.. Was expecting only one file..'
        'Going to use just one'
    fifname = fifs[0]

    # Load data and read event channel
    raw = mne.io.Raw(fpath + fifname, preload=True)
    eves = mne.find_events(raw, stim_channel='STI101', shortest_event=1)

    raw.info['bads'] += ['MEG1421', 'MEG1431', 'MEG2621']
    # Filter the data for ERPs
    raw.filter(l_freq=1.0, h_freq=20, l_trans_bandwidth=0.15,
               picks=np.arange(0, 308, 1))

    # raw.apply_proj()
    fs = raw.info['sfreq']
    # SSP for blinks
    blinks = find_blinks(raw, ch_name='EOG061')
    epochs_blinks = mne.Epochs(raw, blinks, 998, tmin=-0.25,
                               tmax=0.25, proj=True,
                               baseline=(-0.25, 0),
                               reject=dict(grad=8000e-13,
                                           mag=8e-12))
    blink_projs = compute_proj_epochs(epochs_blinks, n_grad=2,
                                      n_mag=2, n_eeg=0,
                                      verbose='DEBUG')
    raw.add_proj(blink_projs)
    evokeds = []
    condnames = ['LL', 'RL', 'LR', 'RR']
    condlists = [[1, 2, 5, 6], [9, 10, 13, 14], [3, 4, 7, 8], [11, 12, 15, 16]]

    for k, condstem in enumerate(condnames):
        condlist = condlists[k]
        print 'Running Subject', subj, 'Condition', condstem

        # Epoching events of type
        epochs = mne.Epochs(raw, eves, condlist, tmin=-0.3, proj=True,
                            tmax=1.0, baseline=(-0.3, 0.0), name=condstem,
                            reject=dict(grad=5000e-13, mag=5e-12))
        evokeds += [epochs.average(), ]

avename = subj + '_ITDadapt-ave.fif'
mne.write_evokeds(fpath + avename, evokeds)
