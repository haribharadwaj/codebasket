import mne
import numpy as np
import os
import fnmatch
from anlffr import spectral
from anlffr.preproc import find_blinks
from mne.preprocessing.ssp import compute_proj_epochs
import pylab as pl
from scipy.io import savemat


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
froot = '/autofs/cluster/transcend/hari/ASSRold/'
saveResults = True
subjlist = ['011302', ]
ch = range(306)  # Channels of interest
mags = range(2, 306, 3)
grads = range(0, 306, 3) + range(1, 306, 3)
freqs = np.arange(5, 500, 2)  # define frequencies of interest
n_cycles = freqs / float(3)  # different number of cycle per frequency
n_cycles[freqs < 15] = 2

ASSR25 = True  # Set false for ASSR43
sss = True
for subj in subjlist:

    fpath = froot + subj + '/'
    if sss:
        ssstag = '_sss'
    else:
        ssstag = ''

    l_freq = 2.0
    if ASSR25:
        condlist = [1, ]
        condstem = 'ASSR25'
    else:
        condlist = [2, ]
        condstem = 'ASSR43'

    print 'Running Subject', subj, 'Condition', condstem

    save_raw_name = subj + '_' + condstem + ssstag + '-epo.fif'

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
        fifs = fnmatch.filter(os.listdir(fpath), subj + '*raw' +
                              ssstag + '.fif')
        print 'No pre-epoched data found, looking for raw files'
        print 'Viola!', len(fifs),  'files found!'
        for k, fif in enumerate(fifs):
            fifs[k] = fpath + fif
        # Load data and read event channel
        raw = mne.io.Raw(fifs, preload=True, add_eeg_ref=False)
        eves = mne.find_events(raw, stim_channel='STI101',
                               shortest_event=1)
        if not sss:
            raw.info['bads'] += ['MEG2021', ]

        # Filter the data for SSRs
        raw.filter(l_freq=l_freq, h_freq=144, l_trans_bandwidth=0.15,
                   picks=np.arange(0, 306, 1))

        # SSP for blinks
        blinks = find_blinks(raw, ch_name='EOG062')
        epochs_blinks = mne.Epochs(raw, blinks, 998, tmin=-0.25,
                                   tmax=0.25, proj=True,
                                   baseline=(-0.25, 0),
                                   reject=dict(grad=8000e-13,
                                               mag=8e-12))
        blink_projs = compute_proj_epochs(epochs_blinks, n_grad=2,
                                          n_mag=2, n_eeg=0,
                                          verbose='DEBUG')
        raw.add_proj(blink_projs)

        # SSP for cardiac
        ekg_name = 'ECG063' if 'ECG063' in raw.ch_names else 'MEG1421'
        ekg_thresh = 200e-6 if 'ECG063' in raw.ch_names else 1.0e-12
        qrs = find_blinks(raw, ch_name=ekg_name, h_freq=100.0,
                          thresh=ekg_thresh, event_id=999)
        epochs_qrs = mne.Epochs(raw, qrs, 999, tmin=-0.1,
                                tmax=0.1, proj=True,
                                baseline=(-0.1, 0),
                                reject=dict(grad=8000e-13,
                                            mag=8e-12))
        qrs_projs = compute_proj_epochs(epochs_qrs, n_grad=2,
                                        n_mag=2, n_eeg=0,
                                        verbose='DEBUG')
        raw.add_proj(qrs_projs)
        useProj = True

        # Epoching events of type
        epochs = mne.Epochs(raw, eves, condlist, tmin=-0.1, proj=useProj,
                            tmax=1.3, baseline=(-0.1, 0.),
                            reject=dict(grad=5000e-13, mag=4e-12))

        x = epochs.get_data()
        if saveResults:
            epochs.save(fpath + save_raw_name)
            evoked = epochs.average()
            save_evoked_name = subj + '_' + condstem + ssstag + '-ave.fif'
            evoked.save(fpath + save_evoked_name)

    # Calculate power, plv
    if not preEpoched:
        Fs = raw.info['sfreq']
        times = epochs.times

    #########################################################################
    # Fourier domain stuff
    pl.figure()
    params = dict(Fs=Fs, fpass=[5, 140], tapers=[1, 1], itc=0)
    for badname in epochs.info['bads']:
        bad_ind = epochs.info['ch_names'].index(badname)
        if bad_ind in mags:
            mags.remove(bad_ind)

    y = x[:, mags, :].transpose((1, 0, 2))
    plv, f = spectral.mtplv(y, params, verbose='DEBUG')
    pl.plot(f, plv.T, linewidth=2)
    pl.xlabel('Frequency (Hz)', fontsize=16)
    pl.ylabel('Intertrial PLV', fontsize=16)
    pl.title('MEG Magnetometers', fontsize=16)
    pl.xlim([15, 90])
    pl.show()
    figname_mag = subj + '_' + condstem + ssstag + '_mag-results.pdf'
    pl.savefig(fpath + figname_mag)
    lout = mne.find_layout(epochs.info)
    pos = lout.pos[mags]
    if ASSR25:
        f_AM = 25.0
    else:
        f_AM = 43.0
    ind_AM = np.argmin(np.abs(f - f_AM))
    pl.figure()
    mne.viz.plot_topomap(plv[:, ind_AM], pos, sensors='ok')
    pl.show()
    figname_mag_topo = subj + '_' + condstem + ssstag + '_mag_topo-results.pdf'
    pl.savefig(fpath + figname_mag_topo)

    res = dict(plv=plv, f=f, f_AM=f_AM, subj=subj, mags=mags)
    save_res_name = subj + '_' + condstem + ssstag + '_mag-results.mat'
    savemat(fpath + save_res_name, mdict=res)

    for badname in epochs.info['bads']:
        bad_ind = epochs.info['ch_names'].index(badname)
        if bad_ind in grads:
            grads.remove(bad_ind)

    pl.figure()
    y = x[:, grads, :].transpose((1, 0, 2))
    plv, f = spectral.mtplv(y, params, verbose='DEBUG')
    pl.plot(f, plv.T, linewidth=2)
    pl.xlabel('Frequency (Hz)', fontsize=16)
    pl.ylabel('Intertrial PLV', fontsize=16)
    pl.title('MEG Gradiometers', fontsize=16)
    pl.xlim([15, 90])
    pl.show()
    figname_grad = subj + '_' + condstem + ssstag + '_grad-results.pdf'
    pl.savefig(fpath + figname_grad)
    lout = mne.find_layout(epochs.info)
    pos = lout.pos[grads]
    if ASSR25:
        f_AM = 25.0
    else:
        f_AM = 43.0
    ind_AM = np.argmin(np.abs(f - f_AM))
    pl.figure()
    mne.viz.plot_topomap(plv[:, ind_AM], pos, sensors='ok')
    pl.show()
    figname_grad_topo = (subj + '_' + condstem + ssstag +
                         '_grad_topo-results.pdf')
    pl.savefig(fpath + figname_grad_topo)

    res = dict(plv=plv, f=f, f_AM=f_AM, subj=subj, grads=grads)
    save_res_name = subj + '_' + condstem + ssstag + '_grad-results.mat'
    savemat(fpath + save_res_name, mdict=res)
