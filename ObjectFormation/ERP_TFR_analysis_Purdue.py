import mne
import numpy as np
import os
import fnmatch
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs
from mne.cov import compute_covariance
from mne.time_frequency import tfr_multitaper

# Adding Files and locations
# froot = '/Users/hari/Documents/Data/ObjectFormation/'
# froot = '/autofs/cluster/transcend/hari/ObjectFormation/'

froot = '/home/hari/Data/ObjectFormation/'

# subjlist = ['095801', '096603']  # Need to use mask for trigger channel
# subjlist = ['054401', ]  # This crashed when finding blinks...

# Done
# subjlist = ['035201', '038301', '038302', '039001', '042201', '092002',
#            '096301', '096302', '053001', '030801', '032901', '032902',
#            '013703', '014002', '063101', '075401', '011302', '010401',
#            '096901', '096902', '097201', '097301', '097601', '097701',
#            '098001', '098002', '098101', '098501']

subjlist = ['011201', '011202', '011302', '048102', '052402', '052901',
            '052902', '082601', '082802', '082901', '085701', '086901',
            '087401', '089401', '089402', '092301', '093101', '093302',
            '093901', '097901']
para = 'object'
epochs = []
sss = True
eog = False
ekg = False
doTFR = True
saveEpochs = False
saveCov = False
saveAve = True
for subj in subjlist:

    fpath = froot + subj + '/'
    visits = fnmatch.filter(os.listdir(fpath), '?')
    if len(visits) > 0:
        visit = visits[0]
        fpath += visit + '/'
    # These are so that the generated files are organized better
    respath = fpath + 'RES/'
    if (not os.path.isdir(respath)):
            os.mkdir(respath)
    if sss:
        ssstag = '_sss'
    else:
        ssstag = ''

    fifs = fnmatch.filter(os.listdir(fpath), subj + '_' + para +
                          '_decim_?_raw' + ssstag + '.fif')
    if len(fifs) == 0:
        print 'No decimated files found.. Ucing undecimated files'
        fifs = fnmatch.filter(os.listdir(fpath), subj + '_' + para +
                              '_?_raw' + ssstag + '.fif')
    print 'Viola!', len(fifs),  'files found!'
    if len(fifs) > 1:
        print 'Warning! Using multitple raw files!'
        raws = []
    for k, fif in enumerate(fifs):
        # Load data and read event channel
        temp = mne.io.Raw(fpath + fif, preload=True)
        raws += [temp, ]

    raw = mne.concatenate_raws(raws)
    eves = mne.find_events(raw, stim_channel='STI101', shortest_event=1)

    if not sss:
        raw.info['bads'] += ['MEG1013', 'MEG1623', 'MEG2342', 'MEG2513',
                             'MEG2542']
    # Filter the data for ERPs
    raw.filter(l_freq=1.0, h_freq=90, picks=np.arange(0, 306, 1))

    # raw.apply_proj()
    fs = raw.info['sfreq']
    removeblinks = True

    if removeblinks:
        # SSP for blinks
        blinks = find_blinks(raw, ch_name=['EOG062', ],
                             l_trans_bandwidth='auto')
        blinkname = (fpath + subj + '_' + para + '_blinks_erp' + ssstag +
                     '-eve.fif')
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
    if ('ECG063' in raw.ch_names) and ekg:
        qrschan = 'ECG063'
        thresh = 20e-6
    else:
        qrschan = 'MEG1711'
        thresh = 1e-12
    qrs = find_blinks(raw, ch_name=[qrschan, ], h_freq=100.0, event_id=999,
                      thresh=thresh, l_trans_bandwidth='auto')
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
    condnames = ['coh07', 'coh14', 'coh20']
    condlists = [3, 2, 1]
    eves2 = np.zeros((eves.shape[0]*2, 3), dtype=np.int)
    fs_int = int(raw.info['sfreq'])
    for k, row in enumerate(eves):
        eves2[2*k, :] = row + np.asarray([fs_int, 0, 0])
        eves2[2*k + 1, :] = row + np.asarray([2*fs_int, 0, 0])

    for k, condstem in enumerate(condnames):
        condlist = condlists[k]
        print 'Running Subject', subj, 'Condition', condstem

        # Epoching events of type
        epochs = mne.Epochs(raw, eves2, event_id=condlist, tmin=-0.4,
                            proj=True, tmax=1.0, baseline=(-0.2, 0.0),
                            reject=dict(grad=5000e-13, mag=5e-12))
        evokeds += [epochs.average(), ]
        fstem = subj + ssstag + '_' + para
        if doTFR:
            freqs = np.arange(5., 70., 1.)
            n_cycles = freqs * 0.2
            power, itc = tfr_multitaper(epochs, freqs, n_cycles,
                                        time_bandwidth=2.0, n_jobs=-1)
            fname_pow = fstem + '_pow_' + condstem + '-tfr.h5'
            fname_itc = fstem + '_itc_' + condstem + '-tfr.h5'
            power.save(respath + fname_pow)
            itc.save(respath + fname_itc)

        if saveEpochs:
            fname_epochs = fstem + '_' + condstem + '-epo.fif'
            epochs.save(respath + fname_epochs)

    # Now save overall onset N100
    epochs = mne.Epochs(raw, eves, event_id=condlists, tmin=-0.4,
                        proj=True, tmax=1.0, baseline=(-0.2, 0.0),
                        reject=dict(grad=5000e-13, mag=5e-12))
    evokeds += [epochs.average(), ]
    if saveEpochs:
            fname_epochs = fstem + '_onset-epo.fif'
            epochs.save(respath + fname_epochs)

    if saveAve:
        avename = subj + ssstag + '_' + para + '_collapse-ave.fif'
        mne.write_evokeds(respath + avename, evokeds)

    if saveCov:
        # Compute covatiance
        cov = compute_covariance(epochs, tmin=-0.2, tmax=0.0)
        covname = subj + ssstag + '_' + para + '_collapse-cov.fif'
        cov.save(respath + covname)
