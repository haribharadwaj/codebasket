import mne
import os
import fnmatch
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs


# Adding Files and locations
# froot = '/home/hari/Documents/PythonCodes/ASSRmartinos/'
froot = '/autofs/cluster/transcend/hari/ASSRnew/'

subjlist = ['092201', ]
para = 'assrnew'
epochs = []
sss = True
decimtag = '_decim'
for subj in subjlist:

    fpath = froot + subj + '/'

    if sss:
        ssstag = '_sss'
    else:
        ssstag = ''

    fifs = fnmatch.filter(os.listdir(fpath),
                          subj + '*' + decimtag + 'raw' + ssstag + '.fif')
    print 'Viola!', len(fifs),  'files found!'
    if len(fifs) > 1:
        print 'Warning! Using multitple raw files!'
    for k, fif in enumerate(fifs):
        fifs[k] = fpath + fif

    # Load data and read event channel
    raw = mne.io.Raw(fifs, preload=True)
    eves = mne.find_events(raw, stim_channel='STI101', shortest_event=1)

    if not sss:
        raw.info['bads'] += ['MEG2033', 'MEG0442', 'MEG2343', 'MEG1643']
    # Filter both MEG and EEG data when available the data for ERPs
    raw.filter(l_freq=0.5, h_freq=144, l_trans_bandwidth=0.15,
               picks=None)

    # raw.apply_proj()
    fs = raw.info['sfreq']
    # SSP for blinks
    blinks = find_blinks(raw, ch_name='EOG062')
    blinkname = fpath + subj + '_' + para + '_blinks_erp' + ssstag + '-eve.fif'
    mne.write_events(blinkname, blinks)
    epochs_blinks = mne.Epochs(raw, blinks, 998, tmin=-0.25,
                               tmax=0.25, proj=True,
                               baseline=(-0.25, 0),
                               reject=dict(grad=8000e-13,
                                           mag=8e-12))
    blink_projs = compute_proj_epochs(epochs_blinks, n_grad=2,
                                      n_mag=2, n_eeg=0,
                                      verbose='DEBUG')
    raw.add_proj(blink_projs)

    # SSP for cardiac artifact
    qrs = find_blinks(raw, ch_name='MEG1421', h_freq=100.0, event_id=999,
                      thresh=1e-12)
    qrsname = fpath + subj + '_' + para + '_qrs_erp' + ssstag + '-eve.fif'
    mne.write_events(qrsname, qrs)
    epochs_qrs = mne.Epochs(raw, qrs, 999, tmin=-0.1,
                            tmax=0.25, proj=True,
                            baseline=(-0.1, 0),
                            reject=dict(grad=8000e-13,
                                        mag=8e-12))
    qrs_projs = compute_proj_epochs(epochs_qrs, n_grad=2,
                                    n_mag=2, n_eeg=0,
                                    verbose='DEBUG')
    raw.add_proj(qrs_projs)

    evokeds = []
    condnames = ['Jump25', 'Jump43', 'NoJump25', 'NoJump43']
    condlists = [[1, 3], [2, 4], [5, 7, 9, 11], [6, 8, 10, 12]]

    for k, condstem in enumerate(condnames):
        condlist = condlists[k]
        print 'Running Subject', subj, 'Condition', condstem

        # Epoching events of type
        epochs = mne.Epochs(raw, eves, condlist, tmin=-0.4, proj=True,
                            tmax=1.3, baseline=(-0.3, 0.0), name=condstem,
                            reject=dict(grad=5000e-13, mag=5e-12))
        evokeds += [epochs.average(), ]

avename = subj + ssstag + '_JumpDynamics-ave.fif'
mne.write_evokeds(fpath + avename, evokeds)
